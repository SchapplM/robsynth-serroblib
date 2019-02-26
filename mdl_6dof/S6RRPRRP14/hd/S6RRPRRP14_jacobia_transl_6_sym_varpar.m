% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP14
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP14_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP14_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:53:24
% EndTime: 2019-02-26 21:53:24
% DurationCPUTime: 0.20s
% Computational Cost: add. (256->62), mult. (642->99), div. (0->0), fcn. (825->10), ass. (0->41)
t37 = sin(qJ(2));
t38 = sin(qJ(1));
t41 = cos(qJ(2));
t42 = cos(qJ(1));
t50 = cos(pkin(6));
t47 = t42 * t50;
t28 = t37 * t47 + t38 * t41;
t35 = sin(qJ(5));
t39 = cos(qJ(5));
t27 = t38 * t37 - t41 * t47;
t36 = sin(qJ(4));
t40 = cos(qJ(4));
t34 = sin(pkin(6));
t56 = t34 * t42;
t46 = -t27 * t36 + t40 * t56;
t65 = t28 * t39 + t35 * t46;
t64 = -t28 * t35 + t39 * t46;
t51 = r_i_i_C(3) + qJ(6);
t62 = pkin(5) + r_i_i_C(1);
t43 = t51 * t35 + t62 * t39 + pkin(4);
t63 = pkin(2) + pkin(9);
t61 = pkin(10) + r_i_i_C(2);
t58 = t34 * t38;
t57 = t34 * t41;
t55 = t35 * t36;
t54 = t35 * t37;
t53 = t36 * t39;
t52 = t37 * t39;
t49 = t34 * (pkin(3) + pkin(8));
t48 = t38 * t50;
t45 = t27 * t40 + t36 * t56;
t44 = pkin(4) * t36 - t61 * t40 + qJ(3);
t30 = -t37 * t48 + t42 * t41;
t29 = t42 * t37 + t41 * t48;
t26 = -t36 * t57 + t50 * t40;
t14 = t29 * t36 + t40 * t58;
t13 = -t29 * t40 + t36 * t58;
t11 = t26 * t35 - t34 * t52;
t2 = t14 * t39 + t30 * t35;
t1 = t14 * t35 - t30 * t39;
t3 = [-t38 * pkin(1) + t46 * pkin(4) - t27 * qJ(3) - t63 * t28 + t42 * t49 + t61 * t45 + t51 * t65 + t62 * t64, t51 * (t29 * t39 + t30 * t55) - t63 * t29 + t62 * (-t29 * t35 + t30 * t53) + t44 * t30, t29, -t43 * t13 + t61 * t14, -t62 * t1 + t51 * t2, t1; t42 * pkin(1) + t14 * pkin(4) + t29 * qJ(3) + t51 * t1 + t61 * t13 + t62 * t2 + t63 * t30 + t38 * t49, t62 * (-t27 * t35 + t28 * t53) + t51 * (t27 * t39 + t28 * t55) - t63 * t27 + t44 * t28, t27, t43 * t45 - t46 * t61, -t51 * t64 + t62 * t65, -t65; 0 (t62 * (t35 * t41 + t36 * t52) + t51 * (t36 * t54 - t39 * t41) + t44 * t37 + t63 * t41) * t34, -t57, t61 * t26 + t43 * (-t50 * t36 - t40 * t57) t51 * (t26 * t39 + t34 * t54) - t62 * t11, t11;];
Ja_transl  = t3;
