% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:03:18
% EndTime: 2019-02-26 22:03:18
% DurationCPUTime: 0.14s
% Computational Cost: add. (220->39), mult. (145->46), div. (0->0), fcn. (158->12), ass. (0->32)
t25 = qJ(2) + qJ(3);
t20 = pkin(10) + t25;
t14 = sin(t20);
t15 = cos(t20);
t24 = pkin(11) + qJ(6);
t18 = sin(t24);
t45 = r_i_i_C(2) * t18;
t51 = r_i_i_C(3) * t15 + t14 * t45;
t16 = cos(pkin(11)) * pkin(5) + pkin(4);
t27 = -pkin(9) - qJ(5);
t50 = (r_i_i_C(3) - t27) * t14 + t15 * t16 + pkin(3) * cos(t25);
t19 = cos(t24);
t46 = r_i_i_C(1) * t19;
t30 = (-t16 - t46) * t14 - t15 * t27 - pkin(3) * sin(t25);
t22 = cos(qJ(2)) * pkin(2);
t48 = pkin(1) + t22 + t50;
t28 = sin(qJ(1));
t42 = t51 * t28;
t29 = cos(qJ(1));
t41 = t51 * t29;
t40 = t18 * t29;
t39 = t19 * t29;
t38 = t28 * t18;
t37 = t28 * t19;
t35 = pkin(5) * sin(pkin(11)) + pkin(7) + pkin(8) + qJ(4);
t32 = -sin(qJ(2)) * pkin(2) + t30;
t31 = (-t45 + t46) * t15 + t50;
t4 = t15 * t39 + t38;
t3 = -t15 * t40 + t37;
t2 = -t15 * t37 + t40;
t1 = t15 * t38 + t39;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t48 * t28 + t35 * t29, t29 * t32 + t41, t29 * t30 + t41, t28, t29 * t14, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t35 * t28 + t48 * t29, t28 * t32 + t42, t28 * t30 + t42, -t29, t28 * t14, -r_i_i_C(1) * t1 + t2 * r_i_i_C(2); 0, t22 + t31, t31, 0, -t15 (-r_i_i_C(1) * t18 - r_i_i_C(2) * t19) * t14;];
Ja_transl  = t5;
