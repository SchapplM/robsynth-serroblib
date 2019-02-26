% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR11
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRR11_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR11_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_jacobia_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:20:25
% EndTime: 2019-02-26 21:20:26
% DurationCPUTime: 0.29s
% Computational Cost: add. (602->71), mult. (1494->119), div. (0->0), fcn. (1963->16), ass. (0->57)
t49 = cos(qJ(1));
t70 = cos(pkin(13));
t72 = cos(pkin(6));
t61 = t72 * t70;
t68 = sin(pkin(13));
t74 = sin(qJ(1));
t31 = -t49 * t61 + t68 * t74;
t43 = sin(pkin(6));
t69 = sin(pkin(7));
t63 = t43 * t69;
t71 = cos(pkin(7));
t82 = t31 * t71 + t49 * t63;
t59 = t72 * t68;
t32 = t49 * t59 + t70 * t74;
t46 = sin(qJ(3));
t75 = cos(qJ(3));
t18 = t32 * t75 - t46 * t82;
t64 = t43 * t71;
t27 = t31 * t69 - t49 * t64;
t45 = sin(qJ(4));
t48 = cos(qJ(4));
t10 = t18 * t48 + t27 * t45;
t81 = -t18 * t45 + t27 * t48;
t54 = t49 * t68 + t61 * t74;
t80 = t54 * t71 - t63 * t74;
t17 = t32 * t46 + t75 * t82;
t42 = qJ(5) + qJ(6);
t40 = sin(t42);
t41 = cos(t42);
t79 = (-t10 * t40 + t17 * t41) * r_i_i_C(1) + (-t10 * t41 - t17 * t40) * r_i_i_C(2);
t33 = t49 * t70 - t59 * t74;
t22 = t33 * t75 - t46 * t80;
t51 = t54 * t69 + t64 * t74;
t14 = t22 * t48 + t45 * t51;
t21 = t33 * t46 + t75 * t80;
t5 = -t14 * t40 + t21 * t41;
t6 = t14 * t41 + t21 * t40;
t78 = r_i_i_C(1) * t5 - r_i_i_C(2) * t6;
t58 = t71 * t70;
t60 = t72 * t69;
t25 = t46 * t60 + (t46 * t58 + t68 * t75) * t43;
t30 = -t63 * t70 + t71 * t72;
t16 = t25 * t48 + t30 * t45;
t24 = -t75 * t60 + (t46 * t68 - t58 * t75) * t43;
t77 = (-t16 * t40 + t24 * t41) * r_i_i_C(1) + (-t16 * t41 - t24 * t40) * r_i_i_C(2);
t76 = r_i_i_C(3) + pkin(12) + pkin(11);
t73 = t49 * t43;
t44 = sin(qJ(5));
t67 = pkin(5) * t44 + pkin(10);
t66 = t74 * t43;
t47 = cos(qJ(5));
t39 = pkin(5) * t47 + pkin(4);
t57 = r_i_i_C(1) * t41 - r_i_i_C(2) * t40 + t39;
t55 = r_i_i_C(1) * t40 + r_i_i_C(2) * t41 + t67;
t53 = -t45 * t76 - t48 * t57 - pkin(3);
t13 = t22 * t45 - t48 * t51;
t1 = [-t74 * pkin(1) - t32 * pkin(2) - pkin(3) * t18 - t27 * pkin(9) + qJ(2) * t73 - t10 * t57 - t55 * t17 + t76 * t81, t66, t21 * t53 + t22 * t55, -t13 * t57 + t14 * t76 (-t14 * t44 + t21 * t47) * pkin(5) + t78, t78; t49 * pkin(1) + t33 * pkin(2) + t22 * pkin(3) + pkin(9) * t51 + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + qJ(2) * t66 + t13 * t76 + t14 * t39 + t21 * t67, -t73, t17 * t53 + t18 * t55, t10 * t76 + t57 * t81 (-t10 * t44 + t17 * t47) * pkin(5) + t79, t79; 0, t72, t24 * t53 + t25 * t55, t76 * t16 + t57 * (-t25 * t45 + t30 * t48) (-t16 * t44 + t24 * t47) * pkin(5) + t77, t77;];
Ja_transl  = t1;
