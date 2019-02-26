% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR10_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR10_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:46
% EndTime: 2019-02-26 22:08:47
% DurationCPUTime: 0.14s
% Computational Cost: add. (165->34), mult. (412->57), div. (0->0), fcn. (524->10), ass. (0->30)
t19 = sin(qJ(3));
t22 = cos(qJ(3));
t16 = sin(pkin(11));
t18 = cos(pkin(11));
t26 = r_i_i_C(1) * t16 + r_i_i_C(2) * t18 + qJ(4);
t30 = pkin(3) + r_i_i_C(3) + qJ(5);
t35 = t26 * t19 + t30 * t22 + pkin(2);
t34 = cos(qJ(1));
t17 = sin(pkin(6));
t21 = sin(qJ(1));
t33 = t17 * t21;
t32 = t17 * t22;
t31 = cos(pkin(6));
t29 = t17 * t34;
t28 = t21 * t31;
t27 = t31 * t34;
t25 = t18 * r_i_i_C(1) - t16 * r_i_i_C(2) + pkin(4) + pkin(9);
t20 = sin(qJ(2));
t23 = cos(qJ(2));
t10 = t20 * t27 + t21 * t23;
t1 = t10 * t19 + t22 * t29;
t2 = t10 * t22 - t19 * t29;
t12 = -t20 * t28 + t34 * t23;
t11 = t34 * t20 + t23 * t28;
t9 = t21 * t20 - t23 * t27;
t8 = t31 * t19 + t20 * t32;
t7 = t17 * t20 * t19 - t31 * t22;
t6 = t12 * t22 + t19 * t33;
t5 = t12 * t19 - t21 * t32;
t3 = [-t21 * pkin(1) - t10 * pkin(2) + pkin(8) * t29 - t26 * t1 - t30 * t2 - t25 * t9, -t11 * t35 + t25 * t12, t26 * t6 - t30 * t5, t5, t6, 0; t34 * pkin(1) + t12 * pkin(2) + pkin(8) * t33 + t25 * t11 + t26 * t5 + t30 * t6, t25 * t10 - t35 * t9, -t30 * t1 + t26 * t2, t1, t2, 0; 0 (t25 * t20 + t35 * t23) * t17, t26 * t8 - t30 * t7, t7, t8, 0;];
Ja_transl  = t3;
