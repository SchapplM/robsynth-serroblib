% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRP3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRP3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:51:33
% EndTime: 2019-02-26 19:51:33
% DurationCPUTime: 0.17s
% Computational Cost: add. (215->37), mult. (347->64), div. (0->0), fcn. (438->11), ass. (0->33)
t43 = pkin(5) + r_i_i_C(1);
t18 = pkin(11) + qJ(4);
t16 = sin(t18);
t17 = cos(t18);
t23 = sin(qJ(5));
t25 = cos(qJ(5));
t30 = -t23 * r_i_i_C(2) + t43 * t25 + pkin(4);
t41 = r_i_i_C(3) + qJ(6) + pkin(9);
t42 = t41 * t16 + t30 * t17 + cos(pkin(11)) * pkin(3) + pkin(2);
t19 = sin(pkin(10));
t20 = sin(pkin(6));
t40 = t19 * t20;
t24 = sin(qJ(2));
t39 = t20 * t24;
t26 = cos(qJ(2));
t38 = t20 * t26;
t37 = cos(pkin(6));
t36 = cos(pkin(10));
t35 = t19 * t37;
t34 = t20 * t36;
t31 = t37 * t36;
t28 = t25 * r_i_i_C(2) + t43 * t23 + pkin(8) + qJ(3);
t10 = -t24 * t35 + t36 * t26;
t9 = t36 * t24 + t26 * t35;
t8 = t19 * t26 + t24 * t31;
t7 = t19 * t24 - t26 * t31;
t6 = t37 * t16 + t17 * t39;
t5 = t16 * t39 - t37 * t17;
t4 = t10 * t17 + t16 * t40;
t3 = t10 * t16 - t17 * t40;
t2 = -t16 * t34 + t8 * t17;
t1 = t8 * t16 + t17 * t34;
t11 = [0, t28 * t10 - t42 * t9, t9, -t30 * t3 + t41 * t4 (-t9 * t23 - t4 * t25) * r_i_i_C(2) + t43 * (-t4 * t23 + t9 * t25) t3; 0, t28 * t8 - t42 * t7, t7, -t30 * t1 + t41 * t2 (-t2 * t25 - t7 * t23) * r_i_i_C(2) + t43 * (-t2 * t23 + t7 * t25) t1; 1 (t28 * t24 + t42 * t26) * t20, -t38, -t30 * t5 + t41 * t6 (t23 * t38 - t6 * t25) * r_i_i_C(2) + t43 * (-t6 * t23 - t25 * t38) t5;];
Ja_transl  = t11;
