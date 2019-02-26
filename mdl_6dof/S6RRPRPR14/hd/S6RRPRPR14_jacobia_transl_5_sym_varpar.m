% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR14_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR14_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:45:25
% EndTime: 2019-02-26 21:45:25
% DurationCPUTime: 0.11s
% Computational Cost: add. (118->34), mult. (284->51), div. (0->0), fcn. (357->8), ass. (0->30)
t36 = pkin(4) - r_i_i_C(2);
t16 = sin(pkin(6));
t20 = sin(qJ(1));
t35 = t16 * t20;
t22 = cos(qJ(2));
t34 = t16 * t22;
t23 = cos(qJ(1));
t33 = t16 * t23;
t19 = sin(qJ(2));
t32 = t20 * t19;
t31 = t20 * t22;
t30 = t23 * t19;
t29 = t23 * t22;
t28 = r_i_i_C(3) + qJ(5);
t27 = pkin(2) + pkin(9) + r_i_i_C(1);
t26 = t16 * (pkin(3) + pkin(8));
t17 = cos(pkin(6));
t10 = -t17 * t29 + t32;
t18 = sin(qJ(4));
t21 = cos(qJ(4));
t25 = -t10 * t18 + t21 * t33;
t3 = t10 * t21 + t18 * t33;
t24 = t36 * t18 - t28 * t21 + qJ(3);
t13 = -t17 * t32 + t29;
t12 = t17 * t31 + t30;
t11 = t17 * t30 + t31;
t8 = t17 * t18 + t21 * t34;
t2 = t12 * t18 + t21 * t35;
t1 = -t12 * t21 + t18 * t35;
t4 = [-t20 * pkin(1) - t10 * qJ(3) - t27 * t11 + t23 * t26 + t36 * t25 + t28 * t3, -t27 * t12 + t24 * t13, t12, -t36 * t1 + t28 * t2, t1, 0; t23 * pkin(1) + t12 * qJ(3) + t28 * t1 + t27 * t13 + t36 * t2 + t20 * t26, -t27 * t10 + t24 * t11, t10, -t28 * t25 + t36 * t3, -t3, 0; 0 (t24 * t19 + t27 * t22) * t16, -t34, t28 * (t17 * t21 - t18 * t34) - t36 * t8, t8, 0;];
Ja_transl  = t4;
