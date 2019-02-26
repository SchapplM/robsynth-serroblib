% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:54:00
% EndTime: 2019-02-26 21:54:01
% DurationCPUTime: 0.16s
% Computational Cost: add. (288->37), mult. (167->45), div. (0->0), fcn. (177->12), ass. (0->34)
t26 = qJ(2) + pkin(11);
t23 = qJ(4) + t26;
t22 = qJ(5) + t23;
t18 = sin(t22);
t19 = cos(t22);
t27 = sin(qJ(6));
t46 = r_i_i_C(2) * t27;
t54 = pkin(10) + r_i_i_C(3);
t55 = t18 * t46 + t19 * t54;
t52 = t19 * pkin(5) + t54 * t18;
t51 = pkin(3) * cos(t26) + cos(qJ(2)) * pkin(2);
t29 = cos(qJ(6));
t47 = r_i_i_C(1) * t29;
t35 = (-pkin(5) - t47) * t18;
t32 = t35 - pkin(4) * sin(t23);
t17 = pkin(4) * cos(t23);
t50 = pkin(1) + t17 + t51 + t52;
t30 = cos(qJ(1));
t43 = t27 * t30;
t28 = sin(qJ(1));
t42 = t28 * t27;
t41 = t28 * t29;
t40 = t29 * t30;
t39 = t55 * t28;
t37 = t55 * t30;
t34 = -pkin(3) * sin(t26) - sin(qJ(2)) * pkin(2) + t32;
t33 = (-t46 + t47) * t19 + t52;
t31 = t17 + t33;
t24 = -pkin(9) - pkin(8) - qJ(3) - pkin(7);
t5 = t19 * t40 + t42;
t4 = -t19 * t43 + t41;
t3 = -t19 * t41 + t43;
t2 = t19 * t42 + t40;
t1 = [t3 * r_i_i_C(1) + t2 * r_i_i_C(2) - t24 * t30 - t50 * t28, t34 * t30 + t37, t28, t32 * t30 + t37, t30 * t35 + t37, r_i_i_C(1) * t4 - r_i_i_C(2) * t5; t5 * r_i_i_C(1) + t4 * r_i_i_C(2) - t28 * t24 + t50 * t30, t34 * t28 + t39, -t30, t32 * t28 + t39, t28 * t35 + t39, -r_i_i_C(1) * t2 + r_i_i_C(2) * t3; 0, t31 + t51, 0, t31, t33 (-r_i_i_C(1) * t27 - r_i_i_C(2) * t29) * t18;];
Ja_transl  = t1;
