% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP9_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP9_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:50:27
% EndTime: 2019-02-26 21:50:28
% DurationCPUTime: 0.19s
% Computational Cost: add. (242->51), mult. (399->82), div. (0->0), fcn. (505->12), ass. (0->38)
t26 = sin(qJ(2));
t27 = sin(qJ(1));
t29 = cos(qJ(2));
t30 = cos(qJ(1));
t38 = cos(pkin(6));
t36 = t30 * t38;
t11 = t27 * t26 - t29 * t36;
t25 = sin(qJ(5));
t28 = cos(qJ(5));
t12 = t26 * t36 + t27 * t29;
t21 = pkin(11) + qJ(4);
t19 = sin(t21);
t20 = cos(t21);
t23 = sin(pkin(6));
t39 = t23 * t30;
t4 = t12 * t20 - t19 * t39;
t48 = -t11 * t28 + t4 * t25;
t47 = -t11 * t25 - t4 * t28;
t18 = cos(pkin(11)) * pkin(3) + pkin(2);
t34 = t28 * r_i_i_C(1) - t25 * r_i_i_C(2) + pkin(4);
t45 = pkin(10) + r_i_i_C(3);
t46 = t45 * t19 + t34 * t20 + t18;
t42 = t23 * t26;
t41 = t23 * t27;
t40 = t23 * t29;
t37 = t27 * t38;
t35 = t23 * (pkin(3) * sin(pkin(11)) + pkin(8));
t24 = -pkin(9) - qJ(3);
t33 = t25 * r_i_i_C(1) + t28 * r_i_i_C(2) - t24;
t32 = -t12 * t19 - t20 * t39;
t14 = -t26 * t37 + t30 * t29;
t13 = t30 * t26 + t29 * t37;
t10 = t38 * t19 + t20 * t42;
t8 = t14 * t20 + t19 * t41;
t7 = t14 * t19 - t20 * t41;
t2 = t13 * t25 + t8 * t28;
t1 = t13 * t28 - t8 * t25;
t3 = [-t27 * pkin(1) - t4 * pkin(4) + t47 * r_i_i_C(1) + t48 * r_i_i_C(2) + t11 * t24 - t12 * t18 + t30 * t35 + t45 * t32, -t13 * t46 + t33 * t14, t13, -t34 * t7 + t45 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t30 * pkin(1) + t8 * pkin(4) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t13 * t24 + t14 * t18 + t27 * t35 + t45 * t7, -t11 * t46 + t33 * t12, t11, t34 * t32 + t45 * t4, -t48 * r_i_i_C(1) + t47 * r_i_i_C(2), 0; 0 (t33 * t26 + t29 * t46) * t23, -t40, t45 * t10 + t34 * (-t19 * t42 + t38 * t20) (-t10 * t25 - t28 * t40) * r_i_i_C(1) + (-t10 * t28 + t25 * t40) * r_i_i_C(2), 0;];
Ja_transl  = t3;
