% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:19:17
% EndTime: 2019-02-26 22:19:18
% DurationCPUTime: 0.21s
% Computational Cost: add. (455->60), mult. (529->90), div. (0->0), fcn. (662->14), ass. (0->44)
t62 = pkin(11) + r_i_i_C(3);
t38 = sin(qJ(6));
t41 = cos(qJ(6));
t66 = t41 * r_i_i_C(1) - t38 * r_i_i_C(2) + pkin(5);
t39 = sin(qJ(2));
t40 = sin(qJ(1));
t42 = cos(qJ(2));
t43 = cos(qJ(1));
t53 = cos(pkin(6));
t50 = t43 * t53;
t21 = t39 * t50 + t40 * t42;
t36 = qJ(3) + pkin(12);
t33 = qJ(5) + t36;
t31 = sin(t33);
t32 = cos(t33);
t37 = sin(pkin(6));
t54 = t37 * t43;
t10 = t21 * t32 - t31 * t54;
t20 = t40 * t39 - t42 * t50;
t65 = t10 * t38 - t20 * t41;
t64 = -t10 * t41 - t20 * t38;
t28 = pkin(4) * cos(t36) + cos(qJ(3)) * pkin(3);
t26 = pkin(2) + t28;
t63 = t62 * t31 + t66 * t32 + t26;
t57 = t37 * t39;
t56 = t37 * t40;
t55 = t37 * t42;
t27 = pkin(4) * sin(t36) + sin(qJ(3)) * pkin(3);
t52 = t37 * (pkin(8) + t27);
t51 = t40 * t53;
t35 = -pkin(10) - qJ(4) - pkin(9);
t48 = t38 * r_i_i_C(1) + t41 * r_i_i_C(2) - t35;
t9 = -t21 * t31 - t32 * t54;
t47 = t62 * t10 + t66 * t9;
t23 = -t39 * t51 + t43 * t42;
t13 = t23 * t31 - t32 * t56;
t14 = t23 * t32 + t31 * t56;
t46 = -t66 * t13 + t62 * t14;
t19 = t53 * t31 + t32 * t57;
t45 = t62 * t19 + t66 * (-t31 * t57 + t53 * t32);
t22 = t43 * t39 + t42 * t51;
t2 = t14 * t41 + t22 * t38;
t1 = -t14 * t38 + t22 * t41;
t3 = [-t40 * pkin(1) - t10 * pkin(5) + t64 * r_i_i_C(1) + t65 * r_i_i_C(2) + t20 * t35 - t21 * t26 + t43 * t52 + t62 * t9, -t22 * t63 + t48 * t23, -t23 * t27 + t28 * t56 + t46, t22, t46, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t43 * pkin(1) + t14 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t62 * t13 - t22 * t35 + t23 * t26 + t40 * t52, -t20 * t63 + t48 * t21, -t21 * t27 - t28 * t54 + t47, t20, t47, -t65 * r_i_i_C(1) + t64 * r_i_i_C(2); 0 (t48 * t39 + t63 * t42) * t37, -t27 * t57 + t53 * t28 + t45, -t55, t45 (-t19 * t38 - t41 * t55) * r_i_i_C(1) + (-t19 * t41 + t38 * t55) * r_i_i_C(2);];
Ja_transl  = t3;
