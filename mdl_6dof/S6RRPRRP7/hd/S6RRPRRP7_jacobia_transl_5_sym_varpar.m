% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRP7
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
% Datum: 2019-02-26 21:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP7_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP7_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_jacobia_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:49:24
% EndTime: 2019-02-26 21:49:25
% DurationCPUTime: 0.12s
% Computational Cost: add. (111->35), mult. (271->53), div. (0->0), fcn. (328->8), ass. (0->30)
t15 = sin(qJ(2));
t18 = cos(qJ(2));
t33 = pkin(2) + pkin(3);
t20 = t15 * qJ(3) + t33 * t18;
t37 = pkin(1) + t20;
t29 = sin(qJ(4));
t30 = cos(qJ(4));
t8 = t15 * t30 - t18 * t29;
t14 = sin(qJ(5));
t17 = cos(qJ(5));
t22 = t17 * r_i_i_C(1) - t14 * r_i_i_C(2) + pkin(4);
t16 = sin(qJ(1));
t3 = t8 * t16;
t31 = -r_i_i_C(3) - pkin(9);
t7 = t15 * t29 + t18 * t30;
t4 = t7 * t16;
t36 = t22 * t3 - t31 * t4;
t19 = cos(qJ(1));
t24 = t19 * t29;
t25 = t19 * t30;
t5 = -t15 * t24 - t18 * t25;
t6 = -t15 * t25 + t18 * t24;
t35 = -t22 * t6 + t31 * t5;
t34 = -t22 * t7 - t31 * t8;
t32 = pkin(7) - pkin(8);
t23 = -t14 * r_i_i_C(1) - t17 * r_i_i_C(2);
t21 = qJ(3) * t18 - t33 * t15;
t2 = -t16 * t14 - t5 * t17;
t1 = t5 * t14 - t16 * t17;
t9 = [-t31 * t3 - t22 * t4 + (t23 + t32) * t19 - t37 * t16, t21 * t19 - t35, t19 * t15, t35, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; -t5 * pkin(4) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t32 * t16 + t37 * t19 - t31 * t6, t21 * t16 - t36, t16 * t15, t36 (-t4 * t14 + t19 * t17) * r_i_i_C(1) + (-t19 * t14 - t4 * t17) * r_i_i_C(2), 0; 0, t20 - t34, -t18, t34, t23 * t8, 0;];
Ja_transl  = t9;
