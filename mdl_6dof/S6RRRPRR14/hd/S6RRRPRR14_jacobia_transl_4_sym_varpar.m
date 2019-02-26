% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR14_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR14_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_jacobia_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:39
% EndTime: 2019-02-26 22:23:40
% DurationCPUTime: 0.12s
% Computational Cost: add. (104->31), mult. (255->53), div. (0->0), fcn. (319->8), ass. (0->27)
t16 = sin(qJ(3));
t19 = cos(qJ(3));
t27 = r_i_i_C(3) + qJ(4);
t32 = pkin(3) - r_i_i_C(2);
t33 = t27 * t16 + t32 * t19 + pkin(2);
t31 = pkin(9) + r_i_i_C(1);
t15 = sin(pkin(6));
t18 = sin(qJ(1));
t30 = t15 * t18;
t29 = t15 * t19;
t21 = cos(qJ(1));
t28 = t15 * t21;
t26 = cos(pkin(6));
t25 = t18 * t26;
t24 = t21 * t26;
t17 = sin(qJ(2));
t20 = cos(qJ(2));
t10 = t17 * t24 + t18 * t20;
t1 = t10 * t16 + t19 * t28;
t23 = -t10 * t19 + t16 * t28;
t12 = -t17 * t25 + t21 * t20;
t11 = t21 * t17 + t20 * t25;
t9 = t18 * t17 - t20 * t24;
t7 = t15 * t17 * t16 - t26 * t19;
t6 = t12 * t19 + t16 * t30;
t5 = t12 * t16 - t18 * t29;
t2 = [-t18 * pkin(1) - t10 * pkin(2) + pkin(8) * t28 - t27 * t1 + t32 * t23 - t31 * t9, -t11 * t33 + t31 * t12, t27 * t6 - t32 * t5, t5, 0, 0; t21 * pkin(1) + t12 * pkin(2) + pkin(8) * t30 + t31 * t11 + t27 * t5 + t32 * t6, t31 * t10 - t33 * t9, -t32 * t1 - t27 * t23, t1, 0, 0; 0 (t31 * t17 + t33 * t20) * t15, t27 * (t26 * t16 + t17 * t29) - t32 * t7, t7, 0, 0;];
Ja_transl  = t2;
