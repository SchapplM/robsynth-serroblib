% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:31:15
% EndTime: 2019-02-26 22:31:16
% DurationCPUTime: 0.16s
% Computational Cost: add. (273->40), mult. (179->48), div. (0->0), fcn. (190->12), ass. (0->35)
t25 = qJ(2) + qJ(3);
t21 = qJ(4) + t25;
t16 = sin(t21);
t17 = cos(t21);
t23 = pkin(11) + qJ(6);
t18 = sin(t23);
t46 = r_i_i_C(2) * t18;
t51 = r_i_i_C(3) * t17 + t16 * t46;
t14 = cos(pkin(11)) * pkin(5) + pkin(4);
t27 = -pkin(10) - qJ(5);
t50 = t17 * t14 + (r_i_i_C(3) - t27) * t16;
t19 = cos(t23);
t47 = r_i_i_C(1) * t19;
t34 = -t17 * t27 + (-t14 - t47) * t16;
t30 = t34 - pkin(3) * sin(t25);
t15 = pkin(3) * cos(t25);
t22 = cos(qJ(2)) * pkin(2);
t49 = pkin(1) + t15 + t22 + t50;
t28 = sin(qJ(1));
t43 = t51 * t28;
t29 = cos(qJ(1));
t42 = t51 * t29;
t41 = t28 * t18;
t40 = t28 * t19;
t39 = t29 * t18;
t38 = t29 * t19;
t36 = pkin(5) * sin(pkin(11)) + pkin(7) + pkin(8) + pkin(9);
t33 = (-t46 + t47) * t17 + t50;
t32 = -sin(qJ(2)) * pkin(2) + t30;
t31 = t15 + t33;
t4 = t17 * t38 + t41;
t3 = -t17 * t39 + t40;
t2 = -t17 * t40 + t39;
t1 = t17 * t41 + t38;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t49 * t28 + t36 * t29, t29 * t32 + t42, t29 * t30 + t42, t29 * t34 + t42, t29 * t16, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t36 * t28 + t49 * t29, t28 * t32 + t43, t28 * t30 + t43, t28 * t34 + t43, t28 * t16, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, t22 + t31, t31, t33, -t17 (-r_i_i_C(1) * t18 - r_i_i_C(2) * t19) * t16;];
Ja_transl  = t5;
