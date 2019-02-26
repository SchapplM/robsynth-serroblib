% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR9
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR9_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR9_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:33:12
% EndTime: 2019-02-26 21:33:12
% DurationCPUTime: 0.17s
% Computational Cost: add. (177->49), mult. (439->79), div. (0->0), fcn. (558->10), ass. (0->35)
t21 = sin(qJ(5));
t25 = cos(qJ(5));
t20 = sin(qJ(6));
t24 = cos(qJ(6));
t31 = t24 * r_i_i_C(1) - t20 * r_i_i_C(2) + pkin(5);
t37 = pkin(2) + qJ(4);
t42 = pkin(10) + r_i_i_C(3);
t43 = t31 * t21 - t42 * t25 + t37;
t19 = sin(pkin(6));
t22 = sin(qJ(2));
t41 = t19 * t22;
t40 = t19 * t25;
t26 = cos(qJ(2));
t39 = t19 * t26;
t27 = cos(qJ(1));
t38 = t19 * t27;
t36 = -pkin(9) + qJ(3);
t35 = cos(pkin(6));
t23 = sin(qJ(1));
t34 = t23 * t35;
t33 = t27 * t35;
t32 = t19 * (pkin(3) + pkin(4) + pkin(8));
t14 = t22 * t33 + t23 * t26;
t7 = -t14 * t21 + t25 * t38;
t30 = t14 * t25 + t21 * t38;
t29 = t20 * r_i_i_C(1) + t24 * r_i_i_C(2) - t36;
t16 = -t22 * t34 + t27 * t26;
t15 = t27 * t22 + t26 * t34;
t13 = t23 * t22 - t26 * t33;
t12 = t21 * t41 + t35 * t25;
t4 = t16 * t21 + t23 * t40;
t3 = t23 * t19 * t21 - t16 * t25;
t2 = -t15 * t20 + t4 * t24;
t1 = -t15 * t24 - t4 * t20;
t5 = [-t23 * pkin(1) + t29 * t13 - t37 * t14 + t27 * t32 + t42 * t30 + t31 * t7, -t15 * t43 - t29 * t16, t15, t16, -t31 * t3 + t42 * t4, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t27 * pkin(1) + t4 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t36 * t15 + t37 * t16 + t23 * t32 + t42 * t3, -t13 * t43 - t29 * t14, t13, t14, t31 * t30 - t42 * t7 (-t13 * t24 + t7 * t20) * r_i_i_C(1) + (t13 * t20 + t7 * t24) * r_i_i_C(2); 0 (-t29 * t22 + t43 * t26) * t19, -t39, t41, t42 * t12 + t31 * (-t35 * t21 + t22 * t40) (-t12 * t20 + t24 * t39) * r_i_i_C(1) + (-t12 * t24 - t20 * t39) * r_i_i_C(2);];
Ja_transl  = t5;
