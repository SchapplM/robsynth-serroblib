% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR11_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR11_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_jacobia_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:39
% EndTime: 2019-02-26 21:06:39
% DurationCPUTime: 0.25s
% Computational Cost: add. (367->59), mult. (1010->99), div. (0->0), fcn. (1330->14), ass. (0->44)
t36 = cos(qJ(1));
t52 = cos(pkin(12));
t53 = cos(pkin(6));
t45 = t53 * t52;
t51 = sin(pkin(12));
t59 = sin(qJ(1));
t22 = -t36 * t45 + t59 * t51;
t29 = sin(pkin(7));
t32 = cos(pkin(7));
t30 = sin(pkin(6));
t55 = t36 * t30;
t63 = -t22 * t32 - t29 * t55;
t44 = t53 * t51;
t23 = t36 * t44 + t59 * t52;
t34 = sin(qJ(3));
t60 = cos(qJ(3));
t10 = t23 * t60 + t63 * t34;
t19 = -t22 * t29 + t32 * t55;
t33 = sin(qJ(4));
t35 = cos(qJ(4));
t1 = t10 * t33 + t19 * t35;
t62 = t10 * t35 - t19 * t33;
t40 = t36 * t51 + t59 * t45;
t49 = t59 * t30;
t61 = -t29 * t49 + t40 * t32;
t28 = sin(pkin(13));
t31 = cos(pkin(13));
t43 = r_i_i_C(1) * t31 - r_i_i_C(2) * t28 + pkin(4);
t54 = r_i_i_C(3) + qJ(5);
t39 = t54 * t33 + t43 * t35 + pkin(3);
t48 = t29 * t53;
t47 = t32 * t52;
t42 = r_i_i_C(1) * t28 + r_i_i_C(2) * t31 + pkin(10);
t41 = -t30 * t52 * t29 + t53 * t32;
t11 = -t23 * t34 + t63 * t60;
t37 = t40 * t29 + t32 * t49;
t24 = t36 * t52 - t59 * t44;
t18 = t34 * t48 + (t34 * t47 + t60 * t51) * t30;
t14 = t24 * t60 - t61 * t34;
t13 = t24 * t34 + t61 * t60;
t7 = t18 * t33 - t41 * t35;
t6 = t14 * t35 + t37 * t33;
t5 = t14 * t33 - t37 * t35;
t2 = [(t11 * t28 - t31 * t62) * r_i_i_C(1) + (t11 * t31 + t28 * t62) * r_i_i_C(2) - t62 * pkin(4) - t10 * pkin(3) + t11 * pkin(10) - t23 * pkin(2) - t59 * pkin(1) + qJ(2) * t55 - t54 * t1 + t19 * pkin(9), t49, -t39 * t13 + t42 * t14, -t43 * t5 + t54 * t6, t5, 0; (t13 * t28 + t31 * t6) * r_i_i_C(1) + (t13 * t31 - t28 * t6) * r_i_i_C(2) + t6 * pkin(4) + t14 * pkin(3) + t13 * pkin(10) + t24 * pkin(2) + t36 * pkin(1) + qJ(2) * t49 + t54 * t5 + t37 * pkin(9), -t55, t42 * t10 + t39 * t11, -t43 * t1 + t54 * t62, t1, 0; 0, t53, t42 * t18 + t39 * (t60 * t48 + (-t51 * t34 + t60 * t47) * t30) t54 * (t18 * t35 + t41 * t33) - t43 * t7, t7, 0;];
Ja_transl  = t2;
