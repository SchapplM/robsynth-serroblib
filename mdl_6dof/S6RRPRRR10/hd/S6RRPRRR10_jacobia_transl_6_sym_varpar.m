% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR10
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR10_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR10_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:59:17
% EndTime: 2019-02-26 21:59:17
% DurationCPUTime: 0.22s
% Computational Cost: add. (371->58), mult. (526->91), div. (0->0), fcn. (663->14), ass. (0->43)
t24 = cos(pkin(12)) * pkin(3) + pkin(2);
t30 = pkin(12) + qJ(4);
t26 = sin(t30);
t27 = cos(t30);
t38 = cos(qJ(5));
t25 = t38 * pkin(5) + pkin(4);
t31 = qJ(5) + qJ(6);
t28 = sin(t31);
t29 = cos(t31);
t45 = r_i_i_C(1) * t29 - r_i_i_C(2) * t28 + t25;
t55 = r_i_i_C(3) + pkin(11) + pkin(10);
t59 = t55 * t26 + t45 * t27 + t24;
t36 = sin(qJ(2));
t37 = sin(qJ(1));
t39 = cos(qJ(2));
t40 = cos(qJ(1));
t50 = cos(pkin(6));
t47 = t40 * t50;
t18 = t36 * t47 + t37 * t39;
t33 = sin(pkin(6));
t51 = t33 * t40;
t10 = t18 * t27 - t26 * t51;
t17 = t37 * t36 - t39 * t47;
t58 = (-t10 * t28 + t17 * t29) * r_i_i_C(1) + (-t10 * t29 - t17 * t28) * r_i_i_C(2);
t48 = t37 * t50;
t20 = -t36 * t48 + t40 * t39;
t53 = t33 * t37;
t14 = t20 * t27 + t26 * t53;
t19 = t40 * t36 + t39 * t48;
t5 = -t14 * t28 + t19 * t29;
t6 = t14 * t29 + t19 * t28;
t57 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t54 = t33 * t36;
t16 = t50 * t26 + t27 * t54;
t52 = t33 * t39;
t56 = (-t16 * t28 - t29 * t52) * r_i_i_C(1) + (-t16 * t29 + t28 * t52) * r_i_i_C(2);
t35 = sin(qJ(5));
t49 = t35 * pkin(5) + pkin(9) + qJ(3);
t46 = t33 * (pkin(3) * sin(pkin(12)) + pkin(8));
t44 = -t18 * t26 - t27 * t51;
t43 = t28 * r_i_i_C(1) + t29 * r_i_i_C(2) + t49;
t13 = t20 * t26 - t27 * t53;
t1 = [-t37 * pkin(1) - t45 * t10 - t43 * t17 - t18 * t24 + t40 * t46 + t55 * t44, -t19 * t59 + t43 * t20, t19, -t45 * t13 + t55 * t14 (-t14 * t35 + t19 * t38) * pkin(5) + t57, t57; t40 * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t55 * t13 + t14 * t25 + t49 * t19 + t20 * t24 + t37 * t46, -t17 * t59 + t43 * t18, t17, t55 * t10 + t45 * t44 (-t10 * t35 + t17 * t38) * pkin(5) + t58, t58; 0 (t43 * t36 + t59 * t39) * t33, -t52, t55 * t16 + t45 * (-t26 * t54 + t50 * t27) (-t16 * t35 - t38 * t52) * pkin(5) + t56, t56;];
Ja_transl  = t1;
