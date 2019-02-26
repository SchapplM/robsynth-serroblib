% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR4
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
% Datum: 2019-02-26 21:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:55:45
% EndTime: 2019-02-26 21:55:46
% DurationCPUTime: 0.31s
% Computational Cost: add. (482->67), mult. (933->105), div. (0->0), fcn. (1220->14), ass. (0->50)
t76 = pkin(11) + r_i_i_C(3);
t46 = sin(qJ(6));
t50 = cos(qJ(6));
t79 = t50 * r_i_i_C(1) - t46 * r_i_i_C(2) + pkin(5);
t43 = sin(pkin(12));
t48 = sin(qJ(2));
t52 = cos(qJ(2));
t66 = cos(pkin(12));
t33 = -t52 * t43 - t48 * t66;
t45 = cos(pkin(6));
t30 = t33 * t45;
t49 = sin(qJ(1));
t53 = cos(qJ(1));
t60 = -t48 * t43 + t52 * t66;
t17 = -t53 * t30 + t49 * t60;
t42 = qJ(4) + qJ(5);
t40 = sin(t42);
t41 = cos(t42);
t44 = sin(pkin(6));
t67 = t53 * t44;
t10 = t17 * t41 - t40 * t67;
t57 = t60 * t45;
t16 = t49 * t33 + t53 * t57;
t78 = t10 * t46 + t16 * t50;
t77 = -t10 * t50 + t16 * t46;
t51 = cos(qJ(4));
t38 = t51 * pkin(4) + pkin(3);
t55 = t76 * t40 + t79 * t41 + t38;
t73 = t52 * pkin(2);
t70 = t45 * t52;
t68 = t49 * t44;
t47 = sin(qJ(4));
t64 = -t45 * t48 * pkin(2) + (t47 * pkin(4) + pkin(8) + qJ(3)) * t44;
t63 = -t49 * t30 - t53 * t60;
t54 = -pkin(10) - pkin(9);
t61 = t46 * r_i_i_C(1) + t50 * r_i_i_C(2) - t54;
t9 = -t17 * t40 - t41 * t67;
t59 = t76 * t10 + t79 * t9;
t13 = -t40 * t63 - t41 * t68;
t14 = t40 * t68 - t41 * t63;
t58 = -t79 * t13 + t76 * t14;
t29 = t33 * t44;
t25 = -t29 * t41 + t45 * t40;
t56 = t76 * t25 + t79 * (t29 * t40 + t45 * t41);
t39 = pkin(1) + t73;
t28 = t60 * t44;
t19 = t53 * t33 - t49 * t57;
t2 = t14 * t50 - t19 * t46;
t1 = -t14 * t46 - t19 * t50;
t3 = [-t10 * pkin(5) + t77 * r_i_i_C(1) + t78 * r_i_i_C(2) - t16 * t54 - t17 * t38 - t49 * t39 + t64 * t53 + t76 * t9 (-t53 * t48 - t49 * t70) * pkin(2) - t61 * t63 + t55 * t19, t68 (t47 * t63 + t51 * t68) * pkin(4) + t58, t58, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t14 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t76 * t13 + t19 * t54 - t38 * t63 + t53 * t39 + t64 * t49 (-t49 * t48 + t53 * t70) * pkin(2) + t61 * t17 + t55 * t16, -t67 (-t17 * t47 - t51 * t67) * pkin(4) + t59, t59, -t78 * r_i_i_C(1) + t77 * r_i_i_C(2); 0, t55 * t28 - t61 * t29 + t44 * t73, t45 (t29 * t47 + t45 * t51) * pkin(4) + t56, t56 (-t25 * t46 - t28 * t50) * r_i_i_C(1) + (-t25 * t50 + t28 * t46) * r_i_i_C(2);];
Ja_transl  = t3;
