% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRP8_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP8_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_jacobia_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:11:54
% EndTime: 2019-02-26 21:11:54
% DurationCPUTime: 0.12s
% Computational Cost: add. (103->33), mult. (125->42), div. (0->0), fcn. (135->8), ass. (0->30)
t15 = qJ(3) + qJ(4);
t13 = sin(t15);
t39 = pkin(9) + r_i_i_C(3);
t41 = t39 * t13;
t14 = cos(t15);
t40 = t39 * t14;
t37 = sin(qJ(3)) * pkin(3);
t36 = cos(qJ(3)) * pkin(3);
t16 = sin(qJ(5));
t35 = r_i_i_C(2) * t16;
t34 = pkin(1) + pkin(8) + pkin(7);
t18 = sin(qJ(1));
t32 = t16 * t18;
t21 = cos(qJ(1));
t31 = t16 * t21;
t19 = cos(qJ(5));
t30 = t18 * t19;
t29 = t19 * t21;
t28 = t14 * t35;
t27 = -r_i_i_C(1) * t19 - pkin(4);
t26 = t18 * t41 + (pkin(4) * t18 + r_i_i_C(1) * t30) * t14;
t25 = t40 + (t27 + t35) * t13;
t24 = t27 * t14 - t41;
t23 = pkin(4) * t13 + qJ(2) + t37 - t40;
t6 = t21 * t28;
t4 = t13 * t29 - t32;
t3 = t13 * t31 + t30;
t2 = t13 * t30 + t31;
t1 = -t13 * t32 + t29;
t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t34 * t18 + t23 * t21, t18 (-t28 + t36) * t18 + t26, -t18 * t28 + t26, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t23 * t18 + t34 * t21, -t21, t6 + (t24 - t36) * t21, t24 * t21 + t6, r_i_i_C(1) * t3 + r_i_i_C(2) * t4, 0; 0, 0, t25 - t37, t25 (-r_i_i_C(1) * t16 - r_i_i_C(2) * t19) * t14, 0;];
Ja_transl  = t5;
