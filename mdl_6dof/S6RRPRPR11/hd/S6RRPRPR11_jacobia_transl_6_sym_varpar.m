% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR11_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR11_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:43:36
% EndTime: 2019-02-26 21:43:36
% DurationCPUTime: 0.11s
% Computational Cost: add. (164->33), mult. (158->45), div. (0->0), fcn. (176->10), ass. (0->28)
t20 = sin(qJ(2));
t22 = cos(qJ(2));
t27 = pkin(2) + r_i_i_C(3) + pkin(9) + qJ(5) + pkin(8);
t26 = t27 * t22;
t19 = qJ(4) + pkin(10);
t10 = pkin(5) * sin(t19) + sin(qJ(4)) * pkin(4);
t28 = qJ(3) + t10;
t36 = -t28 * t20 - pkin(1) - t26;
t11 = pkin(5) * cos(t19) + cos(qJ(4)) * pkin(4);
t34 = pkin(7) + pkin(3) + t11;
t16 = qJ(6) + t19;
t14 = sin(t16);
t15 = cos(t16);
t21 = sin(qJ(1));
t23 = cos(qJ(1));
t29 = t23 * t20;
t5 = -t21 * t14 + t15 * t29;
t6 = t14 * t29 + t21 * t15;
t33 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t30 = t21 * t20;
t7 = t14 * t23 + t15 * t30;
t8 = -t14 * t30 + t15 * t23;
t32 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
t31 = r_i_i_C(1) * t15;
t25 = r_i_i_C(1) * t14 + r_i_i_C(2) * t15 + t28;
t24 = -t27 * t20 + t25 * t22;
t12 = t22 * t14 * r_i_i_C(2);
t1 = [t8 * r_i_i_C(1) - t7 * r_i_i_C(2) + t36 * t21 + t34 * t23, t24 * t23, t29, -t21 * t10 + t11 * t29 + t33, t23 * t22, t33; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t34 * t21 - t36 * t23, t24 * t21, t30, t10 * t23 + t11 * t30 + t32, t21 * t22, t32; 0, t25 * t20 + t26, -t22, t12 + (-t11 - t31) * t22, t20, -t22 * t31 + t12;];
Ja_transl  = t1;
