% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP11
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP11_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP11_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:51:37
% EndTime: 2019-02-26 21:51:38
% DurationCPUTime: 0.11s
% Computational Cost: add. (140->35), mult. (170->48), div. (0->0), fcn. (188->8), ass. (0->26)
t19 = sin(qJ(2));
t21 = cos(qJ(2));
t26 = pkin(2) + r_i_i_C(3) + qJ(6) + pkin(9) + pkin(8);
t25 = t26 * t21;
t18 = qJ(4) + qJ(5);
t14 = sin(t18);
t10 = pkin(5) * t14 + sin(qJ(4)) * pkin(4);
t27 = qJ(3) + t10;
t34 = -t27 * t19 - pkin(1) - t25;
t15 = cos(t18);
t11 = pkin(5) * t15 + cos(qJ(4)) * pkin(4);
t32 = pkin(7) + pkin(3) + t11;
t20 = sin(qJ(1));
t22 = cos(qJ(1));
t28 = t22 * t19;
t5 = -t20 * t14 + t15 * t28;
t6 = t14 * t28 + t20 * t15;
t31 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t29 = t20 * t19;
t7 = t14 * t22 + t15 * t29;
t8 = -t14 * t29 + t15 * t22;
t30 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
t24 = r_i_i_C(1) * t14 + r_i_i_C(2) * t15 + t27;
t23 = -t26 * t19 + t24 * t21;
t12 = t21 * t14 * r_i_i_C(2);
t1 = [t8 * r_i_i_C(1) - t7 * r_i_i_C(2) + t34 * t20 + t32 * t22, t23 * t22, t28, -t20 * t10 + t11 * t28 + t31, t5 * pkin(5) + t31, t22 * t21; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t32 * t20 - t34 * t22, t23 * t20, t29, t22 * t10 + t11 * t29 + t30, t7 * pkin(5) + t30, t20 * t21; 0, t24 * t19 + t25, -t21, t12 + (-r_i_i_C(1) * t15 - t11) * t21, t12 + (-pkin(5) - r_i_i_C(1)) * t21 * t15, t19;];
Ja_transl  = t1;
