% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRPR3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:26:47
% EndTime: 2019-02-26 20:26:47
% DurationCPUTime: 0.10s
% Computational Cost: add. (124->29), mult. (96->36), div. (0->0), fcn. (108->10), ass. (0->23)
t25 = pkin(8) + r_i_i_C(3);
t9 = qJ(4) + pkin(10);
t7 = cos(t9);
t27 = t25 * t7 - sin(qJ(4)) * pkin(4);
t12 = sin(qJ(6));
t14 = cos(qJ(6));
t17 = r_i_i_C(1) * t14 - r_i_i_C(2) * t12 + pkin(5);
t5 = sin(t9);
t26 = t17 * t7 + t25 * t5 + cos(qJ(4)) * pkin(4);
t10 = qJ(1) + pkin(9);
t6 = sin(t10);
t23 = t12 * t6;
t8 = cos(t10);
t22 = t12 * t8;
t20 = t14 * t6;
t19 = t14 * t8;
t18 = pkin(2) + qJ(5) + pkin(7);
t16 = t5 * pkin(5) + qJ(3) - t27;
t4 = t5 * t19 - t23;
t3 = t5 * t22 + t20;
t2 = t5 * t20 + t22;
t1 = -t5 * t23 + t19;
t11 = [-sin(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t18 * t6 + t16 * t8, 0, t6, t26 * t6, t8, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; cos(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t18 * t8 + t16 * t6, 0, -t8, -t26 * t8, t6, r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, 1, 0, -t17 * t5 + t27, 0 (-r_i_i_C(1) * t12 - r_i_i_C(2) * t14) * t7;];
Ja_transl  = t11;
