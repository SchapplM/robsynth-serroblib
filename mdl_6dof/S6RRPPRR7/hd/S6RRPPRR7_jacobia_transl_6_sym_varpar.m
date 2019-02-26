% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR7
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
% Datum: 2019-02-26 21:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:31:56
% EndTime: 2019-02-26 21:31:57
% DurationCPUTime: 0.14s
% Computational Cost: add. (184->50), mult. (448->78), div. (0->0), fcn. (569->10), ass. (0->35)
t40 = pkin(10) + r_i_i_C(3);
t17 = sin(pkin(6));
t20 = sin(qJ(2));
t39 = t17 * t20;
t24 = cos(qJ(2));
t38 = t17 * t24;
t21 = sin(qJ(1));
t37 = t21 * t17;
t25 = cos(qJ(1));
t36 = t25 * t17;
t35 = pkin(4) + qJ(3);
t34 = cos(pkin(6));
t33 = pkin(2) + pkin(3) + pkin(9);
t32 = t17 * (pkin(8) - qJ(4));
t31 = t21 * t34;
t30 = t25 * t34;
t18 = sin(qJ(6));
t22 = cos(qJ(6));
t29 = t22 * r_i_i_C(1) - t18 * r_i_i_C(2) + pkin(5);
t11 = t21 * t20 - t24 * t30;
t19 = sin(qJ(5));
t23 = cos(qJ(5));
t28 = -t11 * t19 + t23 * t36;
t4 = t11 * t23 + t19 * t36;
t27 = t18 * r_i_i_C(1) + t22 * r_i_i_C(2) + t33;
t26 = t40 * t19 + t29 * t23 + t35;
t14 = -t20 * t31 + t25 * t24;
t13 = t25 * t20 + t24 * t31;
t12 = t20 * t30 + t21 * t24;
t10 = t34 * t19 + t23 * t38;
t8 = t13 * t23 - t19 * t37;
t7 = t13 * t19 + t23 * t37;
t2 = t14 * t18 + t8 * t22;
t1 = t14 * t22 - t8 * t18;
t3 = [-t21 * pkin(1) - t35 * t11 - t27 * t12 + t25 * t32 + t40 * t28 - t29 * t4, -t27 * t13 + t26 * t14, t13, -t37, -t29 * t7 + t40 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t25 * pkin(1) + t8 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t35 * t13 + t33 * t14 + t21 * t32 + t40 * t7, -t27 * t11 + t26 * t12, t11, t36, t29 * t28 + t40 * t4 (t12 * t22 - t4 * t18) * r_i_i_C(1) + (-t12 * t18 - t4 * t22) * r_i_i_C(2); 0 (t26 * t20 + t27 * t24) * t17, -t38, -t34, -t40 * t10 + t29 * (t19 * t38 - t34 * t23) (t10 * t18 + t22 * t39) * r_i_i_C(1) + (t10 * t22 - t18 * t39) * r_i_i_C(2);];
Ja_transl  = t3;
