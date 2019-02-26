% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR11
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
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR11_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR11_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_jacobia_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:21:49
% EndTime: 2019-02-26 22:21:50
% DurationCPUTime: 0.12s
% Computational Cost: add. (104->31), mult. (255->53), div. (0->0), fcn. (319->8), ass. (0->27)
t18 = sin(qJ(3));
t21 = cos(qJ(3));
t29 = r_i_i_C(3) + qJ(4);
t34 = pkin(3) + r_i_i_C(1);
t35 = t29 * t18 + t34 * t21 + pkin(2);
t33 = pkin(9) + r_i_i_C(2);
t17 = sin(pkin(6));
t20 = sin(qJ(1));
t32 = t17 * t20;
t31 = t17 * t21;
t23 = cos(qJ(1));
t30 = t17 * t23;
t28 = cos(pkin(6));
t19 = sin(qJ(2));
t22 = cos(qJ(2));
t25 = t23 * t28;
t10 = t19 * t25 + t20 * t22;
t27 = t10 * t21 - t18 * t30;
t26 = t20 * t28;
t1 = t10 * t18 + t21 * t30;
t12 = -t19 * t26 + t23 * t22;
t11 = t23 * t19 + t22 * t26;
t9 = t20 * t19 - t22 * t25;
t7 = t17 * t19 * t18 - t28 * t21;
t6 = t12 * t21 + t18 * t32;
t5 = t12 * t18 - t20 * t31;
t2 = [-t20 * pkin(1) - t10 * pkin(2) + pkin(8) * t30 - t29 * t1 - t34 * t27 - t33 * t9, -t11 * t35 + t33 * t12, t29 * t6 - t34 * t5, t5, 0, 0; t23 * pkin(1) + t12 * pkin(2) + pkin(8) * t32 + t33 * t11 + t29 * t5 + t34 * t6, t33 * t10 - t35 * t9, -t34 * t1 + t29 * t27, t1, 0, 0; 0 (t33 * t19 + t35 * t22) * t17, t29 * (t28 * t18 + t19 * t31) - t34 * t7, t7, 0, 0;];
Ja_transl  = t2;
