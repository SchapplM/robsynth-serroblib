% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPRRPR2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRPR2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:56
% EndTime: 2019-02-26 19:40:57
% DurationCPUTime: 0.18s
% Computational Cost: add. (347->46), mult. (972->85), div. (0->0), fcn. (1282->14), ass. (0->46)
t51 = sin(pkin(12));
t52 = sin(pkin(11));
t41 = t52 * t51;
t55 = cos(pkin(12));
t56 = cos(pkin(11));
t49 = t56 * t55;
t58 = cos(pkin(6));
t33 = -t58 * t49 + t41;
t53 = sin(pkin(7));
t54 = sin(pkin(6));
t46 = t54 * t53;
t57 = cos(pkin(7));
t62 = t33 * t57 + t56 * t46;
t43 = t52 * t55;
t47 = t56 * t51;
t34 = t58 * t43 + t47;
t42 = t52 * t54;
t61 = t34 * t57 - t53 * t42;
t60 = t54 * t55 * t57 + t53 * t58;
t59 = cos(qJ(3));
t50 = -pkin(4) - pkin(10) - r_i_i_C(3);
t48 = t56 * t54;
t45 = t54 * t51;
t23 = sin(qJ(6));
t26 = cos(qJ(6));
t40 = t23 * r_i_i_C(1) + t26 * r_i_i_C(2) + qJ(5);
t39 = t26 * r_i_i_C(1) - t23 * r_i_i_C(2) + pkin(5) + pkin(9);
t24 = sin(qJ(4));
t27 = cos(qJ(4));
t35 = -t40 * t24 + t50 * t27 - pkin(3);
t32 = -t55 * t46 + t58 * t57;
t29 = t34 * t53 + t57 * t42;
t28 = t33 * t53 - t57 * t48;
t25 = sin(qJ(3));
t19 = -t58 * t41 + t49;
t18 = t58 * t47 + t43;
t14 = t60 * t25 + t59 * t45;
t13 = t25 * t45 - t60 * t59;
t9 = t14 * t24 - t32 * t27;
t8 = t19 * t59 - t61 * t25;
t7 = t19 * t25 + t61 * t59;
t6 = t18 * t59 - t62 * t25;
t5 = t18 * t25 + t62 * t59;
t3 = t8 * t24 - t29 * t27;
t1 = t6 * t24 - t28 * t27;
t2 = [0, t42, t35 * t7 + t39 * t8, t40 * (t29 * t24 + t8 * t27) + t50 * t3, t3 (-t7 * t23 + t3 * t26) * r_i_i_C(1) + (-t3 * t23 - t7 * t26) * r_i_i_C(2); 0, -t48, t35 * t5 + t39 * t6, t40 * (t28 * t24 + t6 * t27) + t50 * t1, t1 (t1 * t26 - t5 * t23) * r_i_i_C(1) + (-t1 * t23 - t5 * t26) * r_i_i_C(2); 1, t58, t35 * t13 + t39 * t14, t50 * t9 + t40 * (t14 * t27 + t32 * t24) t9 (-t13 * t23 + t9 * t26) * r_i_i_C(1) + (-t13 * t26 - t9 * t23) * r_i_i_C(2);];
Ja_transl  = t2;
