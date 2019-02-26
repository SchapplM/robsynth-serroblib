% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR9
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
% Datum: 2019-02-26 21:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR9_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR9_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_jacobia_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:05:28
% EndTime: 2019-02-26 21:05:28
% DurationCPUTime: 0.22s
% Computational Cost: add. (249->58), mult. (590->97), div. (0->0), fcn. (767->14), ass. (0->49)
t35 = sin(qJ(4));
t60 = t35 * pkin(4) + pkin(9);
t38 = cos(qJ(4));
t24 = t38 * pkin(4) + pkin(3);
t27 = qJ(4) + pkin(13);
t25 = sin(t27);
t26 = cos(t27);
t59 = t26 * r_i_i_C(1) - t25 * r_i_i_C(2) + t24;
t57 = r_i_i_C(3) + qJ(5) + pkin(10);
t56 = cos(qJ(3));
t29 = sin(pkin(7));
t36 = sin(qJ(3));
t55 = t29 * t36;
t30 = sin(pkin(6));
t31 = cos(pkin(12));
t54 = t30 * t31;
t32 = cos(pkin(7));
t53 = t32 * t36;
t28 = sin(pkin(12));
t37 = sin(qJ(1));
t52 = t37 * t28;
t51 = t37 * t30;
t50 = t37 * t31;
t39 = cos(qJ(1));
t49 = t39 * t28;
t48 = t39 * t30;
t47 = t39 * t31;
t46 = t30 * qJ(2);
t45 = t29 * t56;
t44 = t32 * t56;
t43 = t30 * t45;
t33 = cos(pkin(6));
t17 = -t33 * t47 + t52;
t11 = -t17 * t29 + t32 * t48;
t41 = t33 * t50 + t49;
t40 = t41 * t32;
t18 = t33 * t49 + t50;
t4 = -t17 * t53 + t18 * t56 - t48 * t55;
t13 = t41 * t29 + t32 * t51;
t3 = t17 * t44 + t18 * t36 + t39 * t43;
t19 = -t33 * t52 + t47;
t16 = -t29 * t54 + t33 * t32;
t10 = t33 * t55 + (t56 * t28 + t31 * t53) * t30;
t9 = t30 * t28 * t36 - t33 * t45 - t44 * t54;
t8 = t19 * t56 + (t29 * t51 - t40) * t36;
t7 = t19 * t36 - t37 * t43 + t56 * t40;
t2 = t13 * t25 + t8 * t26;
t1 = t13 * t26 - t8 * t25;
t5 = [-t37 * pkin(1) - t18 * pkin(2) - t57 * t3 + t39 * t46 - t59 * t4 + (t25 * r_i_i_C(1) + t26 * r_i_i_C(2) + t60) * t11, t51, t57 * t8 - t59 * t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2) + (t13 * t38 - t35 * t8) * pkin(4), t7, 0; t39 * pkin(1) + t19 * pkin(2) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t60 * t13 + t8 * t24 + t37 * t46 + t57 * t7, -t48, -t3 * t59 + t57 * t4 (-t11 * t26 - t4 * t25) * r_i_i_C(1) + (t11 * t25 - t4 * t26) * r_i_i_C(2) + (-t11 * t38 - t4 * t35) * pkin(4), t3, 0; 0, t33, t57 * t10 - t59 * t9 (-t10 * t25 + t16 * t26) * r_i_i_C(1) + (-t10 * t26 - t16 * t25) * r_i_i_C(2) + (-t10 * t35 + t16 * t38) * pkin(4), t9, 0;];
Ja_transl  = t5;
