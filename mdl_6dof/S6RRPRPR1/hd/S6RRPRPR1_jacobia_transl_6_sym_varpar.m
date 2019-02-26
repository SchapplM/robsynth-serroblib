% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:37:50
% EndTime: 2019-02-26 21:37:50
% DurationCPUTime: 0.12s
% Computational Cost: add. (214->39), mult. (140->46), div. (0->0), fcn. (153->12), ass. (0->32)
t24 = qJ(2) + pkin(10);
t20 = qJ(4) + t24;
t15 = sin(t20);
t16 = cos(t20);
t23 = pkin(11) + qJ(6);
t18 = sin(t23);
t44 = r_i_i_C(2) * t18;
t48 = r_i_i_C(3) * t16 + t15 * t44;
t17 = cos(pkin(11)) * pkin(5) + pkin(4);
t26 = -pkin(9) - qJ(5);
t47 = t16 * t17 + (r_i_i_C(3) - t26) * t15;
t35 = pkin(3) * cos(t24) + cos(qJ(2)) * pkin(2);
t46 = pkin(1) + t35 + t47;
t19 = cos(t23);
t45 = r_i_i_C(1) * t19;
t27 = sin(qJ(1));
t41 = t48 * t27;
t28 = cos(qJ(1));
t40 = t48 * t28;
t39 = t18 * t28;
t38 = t19 * t28;
t37 = t27 * t18;
t36 = t27 * t19;
t33 = pkin(5) * sin(pkin(11)) + pkin(7) + pkin(8) + qJ(3);
t31 = -t16 * t26 + (-t17 - t45) * t15;
t30 = (-t44 + t45) * t16 + t47;
t29 = -pkin(3) * sin(t24) - sin(qJ(2)) * pkin(2) + t31;
t4 = t16 * t38 + t37;
t3 = -t16 * t39 + t36;
t2 = -t16 * t36 + t39;
t1 = t16 * t37 + t38;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t46 * t27 + t33 * t28, t29 * t28 + t40, t27, t31 * t28 + t40, t28 * t15, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t33 * t27 + t46 * t28, t29 * t27 + t41, -t28, t31 * t27 + t41, t27 * t15, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, t30 + t35, 0, t30, -t16 (-r_i_i_C(1) * t18 - r_i_i_C(2) * t19) * t15;];
Ja_transl  = t5;
