% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR6_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR6_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:50:04
% EndTime: 2019-02-26 22:50:04
% DurationCPUTime: 0.20s
% Computational Cost: add. (314->56), mult. (509->91), div. (0->0), fcn. (637->12), ass. (0->43)
t59 = pkin(11) + r_i_i_C(3);
t32 = sin(qJ(5));
t36 = cos(qJ(5));
t63 = t36 * r_i_i_C(1) - t32 * r_i_i_C(2) + pkin(4);
t34 = sin(qJ(2));
t35 = sin(qJ(1));
t38 = cos(qJ(2));
t39 = cos(qJ(1));
t50 = cos(pkin(6));
t48 = t39 * t50;
t21 = t34 * t48 + t35 * t38;
t30 = qJ(3) + qJ(4);
t28 = sin(t30);
t29 = cos(t30);
t31 = sin(pkin(6));
t51 = t31 * t39;
t10 = t21 * t29 - t28 * t51;
t20 = t35 * t34 - t38 * t48;
t62 = t10 * t32 - t20 * t36;
t61 = -t10 * t36 - t20 * t32;
t37 = cos(qJ(3));
t27 = t37 * pkin(3) + pkin(2);
t60 = t59 * t28 + t63 * t29 + t27;
t54 = t31 * t34;
t53 = t31 * t35;
t52 = t31 * t38;
t49 = t35 * t50;
t33 = sin(qJ(3));
t47 = t31 * (pkin(3) * t33 + pkin(8));
t40 = -pkin(10) - pkin(9);
t45 = t32 * r_i_i_C(1) + t36 * r_i_i_C(2) - t40;
t9 = -t21 * t28 - t29 * t51;
t44 = t59 * t10 + t63 * t9;
t23 = -t34 * t49 + t39 * t38;
t13 = t23 * t28 - t29 * t53;
t14 = t23 * t29 + t28 * t53;
t43 = -t63 * t13 + t59 * t14;
t19 = t50 * t28 + t29 * t54;
t42 = t59 * t19 + t63 * (-t28 * t54 + t50 * t29);
t22 = t39 * t34 + t38 * t49;
t2 = t14 * t36 + t22 * t32;
t1 = -t14 * t32 + t22 * t36;
t3 = [-t35 * pkin(1) - t10 * pkin(4) + t61 * r_i_i_C(1) + t62 * r_i_i_C(2) + t20 * t40 - t21 * t27 + t39 * t47 + t59 * t9, -t22 * t60 + t45 * t23 (-t23 * t33 + t37 * t53) * pkin(3) + t43, t43, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t39 * pkin(1) + t14 * pkin(4) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t59 * t13 - t22 * t40 + t23 * t27 + t35 * t47, -t20 * t60 + t45 * t21 (-t21 * t33 - t37 * t51) * pkin(3) + t44, t44, -t62 * r_i_i_C(1) + t61 * r_i_i_C(2), 0; 0 (t45 * t34 + t60 * t38) * t31 (-t33 * t54 + t50 * t37) * pkin(3) + t42, t42 (-t19 * t32 - t36 * t52) * r_i_i_C(1) + (-t19 * t36 + t32 * t52) * r_i_i_C(2), 0;];
Ja_transl  = t3;
