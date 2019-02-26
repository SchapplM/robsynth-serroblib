% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR5
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
% Datum: 2019-02-26 21:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:56:18
% EndTime: 2019-02-26 21:56:19
% DurationCPUTime: 0.22s
% Computational Cost: add. (438->64), mult. (955->103), div. (0->0), fcn. (1255->14), ass. (0->48)
t40 = sin(pkin(12));
t45 = sin(qJ(2));
t49 = cos(qJ(2));
t61 = cos(pkin(12));
t55 = -t45 * t40 + t49 * t61;
t44 = sin(qJ(4));
t48 = cos(qJ(4));
t47 = cos(qJ(5));
t35 = t47 * pkin(5) + pkin(4);
t39 = qJ(5) + qJ(6);
t37 = sin(t39);
t38 = cos(t39);
t57 = r_i_i_C(1) * t38 - r_i_i_C(2) * t37 + t35;
t66 = r_i_i_C(3) + pkin(11) + pkin(10);
t52 = t66 * t44 + t57 * t48 + pkin(3);
t30 = -t49 * t40 - t45 * t61;
t42 = cos(pkin(6));
t27 = t30 * t42;
t46 = sin(qJ(1));
t50 = cos(qJ(1));
t17 = -t50 * t27 + t46 * t55;
t41 = sin(pkin(6));
t62 = t50 * t41;
t10 = t17 * t48 - t44 * t62;
t53 = t55 * t42;
t16 = t46 * t30 + t50 * t53;
t70 = (-t10 * t37 - t16 * t38) * r_i_i_C(1) + (-t10 * t38 + t16 * t37) * r_i_i_C(2);
t58 = -t46 * t27 - t50 * t55;
t63 = t46 * t41;
t14 = t44 * t63 - t48 * t58;
t19 = t50 * t30 - t46 * t53;
t5 = -t14 * t37 - t19 * t38;
t6 = t14 * t38 - t19 * t37;
t69 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t26 = t30 * t41;
t22 = -t26 * t48 + t42 * t44;
t25 = t55 * t41;
t68 = (-t22 * t37 - t25 * t38) * r_i_i_C(1) + (-t22 * t38 + t25 * t37) * r_i_i_C(2);
t67 = t49 * pkin(2);
t65 = t42 * t49;
t43 = sin(qJ(5));
t60 = -t43 * pkin(5) - pkin(9);
t56 = -t17 * t44 - t48 * t62;
t54 = t37 * r_i_i_C(1) + t38 * r_i_i_C(2) - t60;
t36 = pkin(1) + t67;
t28 = t42 * t45 * pkin(2) + (-pkin(8) - qJ(3)) * t41;
t13 = -t44 * t58 - t48 * t63;
t1 = [-t17 * pkin(3) - t57 * t10 + t54 * t16 - t50 * t28 - t46 * t36 + t66 * t56 (-t50 * t45 - t46 * t65) * pkin(2) - t54 * t58 + t52 * t19, t63, -t57 * t13 + t66 * t14 (-t14 * t43 - t19 * t47) * pkin(5) + t69, t69; -pkin(3) * t58 + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t66 * t13 + t14 * t35 + t60 * t19 - t46 * t28 + t50 * t36 (-t46 * t45 + t50 * t65) * pkin(2) + t54 * t17 + t52 * t16, -t62, t66 * t10 + t57 * t56 (-t10 * t43 - t16 * t47) * pkin(5) + t70, t70; 0, t52 * t25 - t54 * t26 + t41 * t67, t42, t66 * t22 + t57 * (t26 * t44 + t42 * t48) (-t22 * t43 - t25 * t47) * pkin(5) + t68, t68;];
Ja_transl  = t1;
