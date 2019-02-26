% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR3_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR3_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:02:01
% EndTime: 2019-02-26 21:02:01
% DurationCPUTime: 0.12s
% Computational Cost: add. (123->25), mult. (144->37), div. (0->0), fcn. (165->8), ass. (0->21)
t13 = cos(qJ(3));
t11 = sin(qJ(3));
t20 = pkin(8) + r_i_i_C(2);
t16 = t20 * t11;
t23 = pkin(3) * t13 + pkin(2) + t16;
t10 = sin(qJ(4));
t12 = cos(qJ(4));
t17 = r_i_i_C(3) + qJ(5);
t21 = pkin(4) + r_i_i_C(1);
t22 = t17 * t10 + t21 * t12 + pkin(3);
t19 = t10 * t13;
t18 = t12 * t13;
t14 = -t22 * t11 + t20 * t13;
t9 = qJ(1) + pkin(10);
t8 = cos(t9);
t7 = sin(t9);
t4 = t7 * t10 + t8 * t18;
t3 = -t7 * t12 + t8 * t19;
t2 = -t8 * t10 + t7 * t18;
t1 = t8 * t12 + t7 * t19;
t5 = [-sin(qJ(1)) * pkin(1) + t8 * pkin(7) - t21 * t2 - t17 * t1 - t23 * t7, 0, t14 * t8, t17 * t4 - t21 * t3, t3, 0; cos(qJ(1)) * pkin(1) + t7 * pkin(7) + t21 * t4 + t17 * t3 + t23 * t8, 0, t14 * t7, -t21 * t1 + t17 * t2, t1, 0; 0, 1, t22 * t13 + t16 (-t21 * t10 + t17 * t12) * t11, t11 * t10, 0;];
Ja_transl  = t5;
