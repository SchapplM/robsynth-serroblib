% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR9_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR9_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:42:23
% EndTime: 2019-02-26 21:42:23
% DurationCPUTime: 0.20s
% Computational Cost: add. (316->52), mult. (445->82), div. (0->0), fcn. (563->14), ass. (0->39)
t20 = cos(pkin(11)) * pkin(3) + pkin(2);
t26 = pkin(11) + qJ(4);
t22 = sin(t26);
t24 = cos(t26);
t19 = cos(pkin(12)) * pkin(5) + pkin(4);
t25 = pkin(12) + qJ(6);
t21 = sin(t25);
t23 = cos(t25);
t38 = t23 * r_i_i_C(1) - t21 * r_i_i_C(2) + t19;
t48 = r_i_i_C(3) + pkin(10) + qJ(5);
t49 = t48 * t22 + t38 * t24 + t20;
t29 = sin(pkin(6));
t32 = sin(qJ(2));
t47 = t29 * t32;
t33 = sin(qJ(1));
t46 = t29 * t33;
t34 = cos(qJ(2));
t45 = t29 * t34;
t35 = cos(qJ(1));
t44 = t29 * t35;
t43 = cos(pkin(6));
t42 = sin(pkin(12)) * pkin(5) + pkin(9) + qJ(3);
t40 = t35 * t43;
t12 = t32 * t40 + t33 * t34;
t4 = t12 * t24 - t22 * t44;
t41 = t33 * t43;
t39 = t29 * (pkin(3) * sin(pkin(11)) + pkin(8));
t3 = t12 * t22 + t24 * t44;
t37 = t21 * r_i_i_C(1) + t23 * r_i_i_C(2) + t42;
t14 = -t32 * t41 + t35 * t34;
t13 = t35 * t32 + t34 * t41;
t11 = t33 * t32 - t34 * t40;
t10 = t43 * t22 + t24 * t47;
t9 = t22 * t47 - t43 * t24;
t8 = t14 * t24 + t22 * t46;
t7 = t14 * t22 - t24 * t46;
t2 = t13 * t21 + t8 * t23;
t1 = t13 * t23 - t8 * t21;
t5 = [-t33 * pkin(1) - t37 * t11 - t12 * t20 - t48 * t3 + t35 * t39 - t38 * t4, -t13 * t49 + t37 * t14, t13, -t38 * t7 + t48 * t8, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t35 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t42 * t13 + t14 * t20 + t8 * t19 + t33 * t39 + t48 * t7, -t11 * t49 + t37 * t12, t11, -t38 * t3 + t48 * t4, t3 (t11 * t23 - t4 * t21) * r_i_i_C(1) + (-t11 * t21 - t4 * t23) * r_i_i_C(2); 0 (t37 * t32 + t34 * t49) * t29, -t45, t48 * t10 - t38 * t9, t9 (-t10 * t21 - t23 * t45) * r_i_i_C(1) + (-t10 * t23 + t21 * t45) * r_i_i_C(2);];
Ja_transl  = t5;
