% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR4_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR4_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:05:08
% EndTime: 2019-02-26 22:05:08
% DurationCPUTime: 0.11s
% Computational Cost: add. (132->31), mult. (165->45), div. (0->0), fcn. (189->8), ass. (0->26)
t16 = cos(qJ(2));
t13 = sin(qJ(2));
t25 = r_i_i_C(2) + qJ(4) + pkin(8);
t19 = t25 * t13;
t15 = cos(qJ(3));
t7 = pkin(3) * t15 + pkin(2);
t29 = t16 * t7 + pkin(1) + t19;
t22 = r_i_i_C(3) + qJ(5);
t27 = pkin(4) + r_i_i_C(1);
t10 = qJ(3) + pkin(10);
t8 = sin(t10);
t9 = cos(t10);
t28 = t22 * t8 + t27 * t9 + t7;
t12 = sin(qJ(3));
t26 = pkin(3) * t12;
t14 = sin(qJ(1));
t24 = t14 * t16;
t17 = cos(qJ(1));
t23 = t16 * t17;
t20 = pkin(7) + t26;
t18 = -t28 * t13 + t25 * t16;
t4 = t14 * t8 + t9 * t23;
t3 = -t14 * t9 + t8 * t23;
t2 = -t17 * t8 + t9 * t24;
t1 = t17 * t9 + t8 * t24;
t5 = [-t22 * t1 - t29 * t14 + t20 * t17 - t27 * t2, t18 * t17, t22 * t4 - t27 * t3 + (-t12 * t23 + t14 * t15) * pkin(3), t17 * t13, t3, 0; t20 * t14 + t29 * t17 + t22 * t3 + t27 * t4, t18 * t14, t22 * t2 - t27 * t1 + (-t12 * t24 - t15 * t17) * pkin(3), t14 * t13, t1, 0; 0, t28 * t16 + t19 (t22 * t9 - t27 * t8 - t26) * t13, -t16, t13 * t8, 0;];
Ja_transl  = t5;
