% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:48:18
% EndTime: 2019-02-26 22:48:18
% DurationCPUTime: 0.17s
% Computational Cost: add. (303->46), mult. (212->59), div. (0->0), fcn. (228->12), ass. (0->42)
t31 = qJ(4) + qJ(5);
t25 = cos(t31);
t15 = pkin(5) * t25 + cos(qJ(4)) * pkin(4);
t12 = pkin(3) + t15;
t32 = qJ(2) + qJ(3);
t24 = sin(t32);
t26 = cos(t32);
t30 = -pkin(11) - pkin(10) - pkin(9);
t59 = t26 * t12 + (r_i_i_C(3) - t30) * t24;
t29 = cos(qJ(2)) * pkin(2);
t58 = pkin(1) + t29 + t59;
t27 = qJ(6) + t31;
t21 = cos(t27);
t35 = cos(qJ(1));
t46 = t35 * t21;
t20 = sin(t27);
t34 = sin(qJ(1));
t49 = t34 * t20;
t5 = t26 * t49 + t46;
t47 = t35 * t20;
t48 = t34 * t21;
t6 = -t26 * t48 + t47;
t57 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = -t26 * t47 + t48;
t8 = t26 * t46 + t49;
t56 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t23 = sin(t31);
t55 = pkin(5) * t23;
t54 = r_i_i_C(1) * t21;
t53 = r_i_i_C(2) * t20;
t51 = t26 * t34;
t50 = t26 * t35;
t42 = t24 * t53;
t45 = r_i_i_C(3) * t51 + t34 * t42;
t44 = r_i_i_C(3) * t50 + t35 * t42;
t14 = t55 + sin(qJ(4)) * pkin(4);
t43 = t14 + pkin(8) + pkin(7);
t40 = -r_i_i_C(1) * t20 - r_i_i_C(2) * t21;
t39 = -t26 * t30 + (-t12 - t54) * t24;
t38 = (-t53 + t54) * t26 + t59;
t37 = -sin(qJ(2)) * pkin(2) + t39;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t58 * t34 + t43 * t35, t37 * t35 + t44, t39 * t35 + t44, -t14 * t50 + t34 * t15 + t56 (-t23 * t50 + t25 * t34) * pkin(5) + t56, t56; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t43 * t34 + t58 * t35, t37 * t34 + t45, t39 * t34 + t45, -t14 * t51 - t35 * t15 + t57 (-t23 * t51 - t25 * t35) * pkin(5) + t57, t57; 0, t29 + t38, t38 (-t14 + t40) * t24 (t40 - t55) * t24, t40 * t24;];
Ja_transl  = t1;
