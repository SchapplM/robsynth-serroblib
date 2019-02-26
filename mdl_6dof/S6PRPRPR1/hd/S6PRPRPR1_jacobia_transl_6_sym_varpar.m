% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRPR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:46:27
% EndTime: 2019-02-26 19:46:27
% DurationCPUTime: 0.20s
% Computational Cost: add. (281->52), mult. (574->85), div. (0->0), fcn. (753->14), ass. (0->36)
t26 = sin(pkin(11));
t34 = sin(qJ(2));
t37 = cos(qJ(2));
t45 = cos(pkin(11));
t40 = -t34 * t26 + t37 * t45;
t25 = qJ(4) + pkin(12);
t23 = sin(t25);
t24 = cos(t25);
t36 = cos(qJ(4));
t32 = sin(qJ(6));
t35 = cos(qJ(6));
t42 = t35 * r_i_i_C(1) - t32 * r_i_i_C(2) + pkin(5);
t50 = pkin(9) + r_i_i_C(3);
t38 = t36 * pkin(4) + t50 * t23 + t42 * t24 + pkin(3);
t27 = sin(pkin(10));
t28 = sin(pkin(6));
t49 = t27 * t28;
t29 = cos(pkin(10));
t48 = t29 * t28;
t30 = cos(pkin(6));
t47 = t30 * t37;
t19 = -t37 * t26 - t34 * t45;
t17 = t19 * t30;
t7 = -t29 * t17 + t27 * t40;
t43 = -t27 * t17 - t29 * t40;
t41 = t32 * r_i_i_C(1) + t35 * r_i_i_C(2) + pkin(8) + qJ(5);
t39 = t40 * t30;
t33 = sin(qJ(4));
t16 = t19 * t28;
t15 = t40 * t28;
t12 = -t16 * t24 + t30 * t23;
t9 = t29 * t19 - t27 * t39;
t6 = t27 * t19 + t29 * t39;
t4 = t23 * t49 - t24 * t43;
t2 = -t23 * t48 + t7 * t24;
t1 = [0 (-t27 * t47 - t29 * t34) * pkin(2) - t41 * t43 + t38 * t9, t49, t50 * t4 + (t33 * t43 + t36 * t49) * pkin(4) + t42 * (t23 * t43 + t24 * t49) -t9 (-t4 * t32 - t9 * t35) * r_i_i_C(1) + (t9 * t32 - t4 * t35) * r_i_i_C(2); 0 (-t27 * t34 + t29 * t47) * pkin(2) + t41 * t7 + t38 * t6, -t48, t50 * t2 + (-t33 * t7 - t36 * t48) * pkin(4) + t42 * (-t7 * t23 - t24 * t48) -t6 (-t2 * t32 - t6 * t35) * r_i_i_C(1) + (-t2 * t35 + t6 * t32) * r_i_i_C(2); 1, t28 * t37 * pkin(2) + t15 * t38 - t41 * t16, t30, t50 * t12 + (t16 * t33 + t30 * t36) * pkin(4) + t42 * (t16 * t23 + t30 * t24) -t15 (-t12 * t32 - t15 * t35) * r_i_i_C(1) + (-t12 * t35 + t15 * t32) * r_i_i_C(2);];
Ja_transl  = t1;
