% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR13_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR13_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobia_transl_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:03
% EndTime: 2019-02-26 22:23:04
% DurationCPUTime: 0.18s
% Computational Cost: add. (213->69), mult. (584->118), div. (0->0), fcn. (752->12), ass. (0->45)
t30 = sin(pkin(7));
t59 = t30 * pkin(10);
t35 = sin(qJ(2));
t36 = sin(qJ(1));
t38 = cos(qJ(2));
t39 = cos(qJ(1));
t46 = cos(pkin(6));
t42 = t39 * t46;
t20 = t35 * t42 + t36 * t38;
t34 = sin(qJ(3));
t58 = t20 * t34;
t29 = sin(pkin(13));
t57 = t29 * t30;
t32 = cos(pkin(13));
t56 = t30 * t32;
t31 = sin(pkin(6));
t55 = t31 * t36;
t54 = t31 * t39;
t33 = cos(pkin(7));
t53 = t33 * t34;
t37 = cos(qJ(3));
t52 = t33 * t37;
t51 = t34 * t35;
t50 = t34 * t38;
t49 = t35 * t37;
t48 = t37 * t38;
t47 = r_i_i_C(3) + qJ(4);
t45 = t30 * t55;
t44 = t30 * t54;
t43 = t36 * t46;
t41 = t46 * t30;
t40 = t32 * r_i_i_C(1) - t29 * r_i_i_C(2) + pkin(3);
t19 = t36 * t35 - t38 * t42;
t13 = -t19 * t30 + t33 * t54;
t21 = -t39 * t35 - t38 * t43;
t14 = -t21 * t30 + t33 * t55;
t4 = t19 * t53 - t20 * t37 + t34 * t44;
t22 = -t35 * t43 + t39 * t38;
t11 = -t37 * t41 + (-t33 * t48 + t51) * t31;
t10 = t21 * t37 - t22 * t53;
t8 = -t19 * t37 - t20 * t53;
t6 = t22 * t37 + (t21 * t33 + t45) * t34;
t5 = -t21 * t52 + t22 * t34 - t37 * t45;
t1 = t19 * t52 + t37 * t44 + t58;
t2 = [(t13 * t29 + t4 * t32) * r_i_i_C(1) + (t13 * t32 - t4 * t29) * r_i_i_C(2) + t4 * pkin(3) - t20 * pkin(2) - t36 * pkin(1) + pkin(9) * t54 + t47 * (-t58 + (-t19 * t33 - t44) * t37) + t13 * pkin(10) (t10 * t32 + t22 * t57) * r_i_i_C(1) + (-t10 * t29 + t22 * t56) * r_i_i_C(2) + t10 * pkin(3) + t21 * pkin(2) + t22 * t59 + t47 * (t21 * t34 + t22 * t52) -t40 * t5 + t47 * t6, t5, 0, 0; (t14 * t29 + t6 * t32) * r_i_i_C(1) + (t14 * t32 - t6 * t29) * r_i_i_C(2) + t6 * pkin(3) + t22 * pkin(2) + t39 * pkin(1) + pkin(9) * t55 + t47 * t5 + t14 * pkin(10) (t20 * t57 + t8 * t32) * r_i_i_C(1) + (t20 * t56 - t8 * t29) * r_i_i_C(2) + t8 * pkin(3) - t19 * pkin(2) + t20 * t59 + t47 * (-t19 * t34 + t20 * t52) -t40 * t1 - t47 * t4, t1, 0, 0; 0 (t40 * (-t33 * t51 + t48) + t38 * pkin(2) + (r_i_i_C(1) * t29 + r_i_i_C(2) * t32 + pkin(10)) * t35 * t30 + t47 * (t33 * t49 + t50)) * t31, t47 * (t34 * t41 + (t33 * t50 + t49) * t31) - t40 * t11, t11, 0, 0;];
Ja_transl  = t2;
