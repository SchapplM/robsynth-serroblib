% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PPRRRR3
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPRRRR3_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRR3_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobia_transl_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:05
% EndTime: 2019-02-26 19:44:05
% DurationCPUTime: 0.15s
% Computational Cost: add. (154->47), mult. (448->94), div. (0->0), fcn. (589->14), ass. (0->40)
t16 = sin(pkin(13));
t19 = sin(pkin(6));
t43 = t16 * t19;
t24 = cos(pkin(6));
t42 = t16 * t24;
t18 = sin(pkin(7));
t41 = t18 * t19;
t40 = t18 * t24;
t21 = cos(pkin(13));
t39 = t21 * t19;
t38 = t21 * t24;
t22 = cos(pkin(8));
t25 = sin(qJ(4));
t37 = t22 * t25;
t27 = cos(qJ(4));
t36 = t22 * t27;
t23 = cos(pkin(7));
t26 = sin(qJ(3));
t35 = t23 * t26;
t34 = t18 * t39;
t17 = sin(pkin(8));
t33 = t17 * (pkin(10) + r_i_i_C(3));
t15 = sin(pkin(14));
t20 = cos(pkin(14));
t10 = -t16 * t15 + t20 * t38;
t11 = t15 * t38 + t16 * t20;
t28 = cos(qJ(3));
t1 = -t11 * t26 + (t10 * t23 - t34) * t28;
t32 = t1 * t22 + t17 * (-t10 * t18 - t23 * t39);
t12 = -t21 * t15 - t20 * t42;
t13 = -t15 * t42 + t21 * t20;
t29 = t12 * t23 + t16 * t41;
t3 = -t13 * t26 + t29 * t28;
t31 = t17 * (-t12 * t18 + t23 * t43) + t22 * t3;
t5 = t28 * t40 + (t20 * t23 * t28 - t15 * t26) * t19;
t30 = t17 * (-t20 * t41 + t24 * t23) + t22 * t5;
t6 = t26 * t40 + (t15 * t28 + t20 * t35) * t19;
t4 = t13 * t28 + t29 * t26;
t2 = t10 * t35 + t11 * t28 - t26 * t34;
t7 = [0, t43 (t3 * t27 - t4 * t37) * r_i_i_C(1) + (-t3 * t25 - t4 * t36) * r_i_i_C(2) + t3 * pkin(3) + t4 * t33 (-t4 * t25 + t31 * t27) * r_i_i_C(1) + (-t31 * t25 - t4 * t27) * r_i_i_C(2), 0, 0; 0, -t39 (t1 * t27 - t2 * t37) * r_i_i_C(1) + (-t1 * t25 - t2 * t36) * r_i_i_C(2) + t1 * pkin(3) + t2 * t33 (-t2 * t25 + t32 * t27) * r_i_i_C(1) + (-t2 * t27 - t32 * t25) * r_i_i_C(2), 0, 0; 1, t24 (t5 * t27 - t6 * t37) * r_i_i_C(1) + (-t5 * t25 - t6 * t36) * r_i_i_C(2) + t5 * pkin(3) + t6 * t33 (-t6 * t25 + t30 * t27) * r_i_i_C(1) + (-t30 * t25 - t6 * t27) * r_i_i_C(2), 0, 0;];
Ja_transl  = t7;
