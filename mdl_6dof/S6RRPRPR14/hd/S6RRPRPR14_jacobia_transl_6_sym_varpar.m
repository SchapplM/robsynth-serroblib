% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR14_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR14_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:45:25
% EndTime: 2019-02-26 21:45:26
% DurationCPUTime: 0.14s
% Computational Cost: add. (197->48), mult. (488->78), div. (0->0), fcn. (621->10), ass. (0->34)
t19 = sin(pkin(6));
t22 = sin(qJ(2));
t41 = t19 * t22;
t23 = sin(qJ(1));
t40 = t19 * t23;
t26 = cos(qJ(2));
t39 = t19 * t26;
t27 = cos(qJ(1));
t38 = t19 * t27;
t37 = cos(pkin(6));
t36 = pkin(2) + pkin(5) + pkin(9);
t35 = pkin(4) + pkin(10) + r_i_i_C(3);
t34 = t19 * (pkin(3) + pkin(8));
t33 = t23 * t37;
t32 = t27 * t37;
t20 = sin(qJ(6));
t24 = cos(qJ(6));
t31 = t20 * r_i_i_C(1) + t24 * r_i_i_C(2) + qJ(5);
t13 = t23 * t22 - t26 * t32;
t21 = sin(qJ(4));
t25 = cos(qJ(4));
t30 = -t13 * t21 + t25 * t38;
t5 = t13 * t25 + t21 * t38;
t29 = t24 * r_i_i_C(1) - t20 * r_i_i_C(2) + t36;
t28 = t35 * t21 - t31 * t25 + qJ(3);
t16 = -t22 * t33 + t27 * t26;
t15 = t27 * t22 + t26 * t33;
t14 = t22 * t32 + t23 * t26;
t11 = t37 * t21 + t25 * t39;
t4 = t15 * t21 + t25 * t40;
t3 = -t15 * t25 + t21 * t40;
t2 = t16 * t24 + t3 * t20;
t1 = -t16 * t20 + t3 * t24;
t6 = [-t23 * pkin(1) - t13 * qJ(3) - t29 * t14 + t27 * t34 + t35 * t30 + t31 * t5, -t29 * t15 + t28 * t16, t15, -t35 * t3 + t31 * t4, t3, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t27 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t15 * qJ(3) + t3 * qJ(5) + t36 * t16 + t23 * t34 + t35 * t4, -t29 * t13 + t28 * t14, t13, -t31 * t30 + t35 * t5, -t5 (-t14 * t20 - t5 * t24) * r_i_i_C(1) + (-t14 * t24 + t5 * t20) * r_i_i_C(2); 0 (t28 * t22 + t29 * t26) * t19, -t39, t31 * (-t21 * t39 + t37 * t25) - t35 * t11, t11 (t11 * t24 - t20 * t41) * r_i_i_C(1) + (-t11 * t20 - t24 * t41) * r_i_i_C(2);];
Ja_transl  = t6;
