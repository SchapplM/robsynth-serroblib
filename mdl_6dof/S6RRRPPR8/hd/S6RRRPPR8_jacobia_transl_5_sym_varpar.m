% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR8_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR8_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:07:39
% EndTime: 2019-02-26 22:07:39
% DurationCPUTime: 0.12s
% Computational Cost: add. (134->33), mult. (325->54), div. (0->0), fcn. (409->8), ass. (0->27)
t18 = sin(qJ(3));
t21 = cos(qJ(3));
t29 = pkin(3) + pkin(4) - r_i_i_C(2);
t31 = r_i_i_C(1) + qJ(4);
t35 = t31 * t18 + t29 * t21 + pkin(2);
t17 = sin(pkin(6));
t20 = sin(qJ(1));
t34 = t17 * t20;
t33 = t17 * t21;
t23 = cos(qJ(1));
t32 = t17 * t23;
t30 = cos(pkin(6));
t28 = pkin(9) - r_i_i_C(3) - qJ(5);
t27 = t20 * t30;
t26 = t23 * t30;
t19 = sin(qJ(2));
t22 = cos(qJ(2));
t10 = t19 * t26 + t20 * t22;
t1 = t10 * t18 + t21 * t32;
t25 = -t10 * t21 + t18 * t32;
t12 = -t19 * t27 + t23 * t22;
t11 = t23 * t19 + t22 * t27;
t9 = t20 * t19 - t22 * t26;
t7 = t17 * t19 * t18 - t30 * t21;
t6 = t12 * t21 + t18 * t34;
t5 = t12 * t18 - t20 * t33;
t2 = [-t20 * pkin(1) - t10 * pkin(2) + pkin(8) * t32 - t31 * t1 + t29 * t25 - t28 * t9, -t11 * t35 + t28 * t12, -t29 * t5 + t31 * t6, t5, -t11, 0; t23 * pkin(1) + t12 * pkin(2) + pkin(8) * t34 + t28 * t11 + t29 * t6 + t31 * t5, t28 * t10 - t35 * t9, -t29 * t1 - t31 * t25, t1, -t9, 0; 0 (t28 * t19 + t35 * t22) * t17, t31 * (t30 * t18 + t19 * t33) - t29 * t7, t7, t17 * t22, 0;];
Ja_transl  = t2;
