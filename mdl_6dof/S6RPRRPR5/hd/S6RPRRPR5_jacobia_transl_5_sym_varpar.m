% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR5
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
% Datum: 2019-02-26 21:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR5_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR5_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:03:12
% EndTime: 2019-02-26 21:03:12
% DurationCPUTime: 0.09s
% Computational Cost: add. (118->20), mult. (73->21), div. (0->0), fcn. (78->7), ass. (0->18)
t34 = r_i_i_C(3) + qJ(5);
t17 = pkin(10) + qJ(3);
t15 = qJ(4) + t17;
t12 = sin(t15);
t13 = cos(t15);
t20 = (pkin(4) - r_i_i_C(2)) * t13 + t34 * t12;
t33 = pkin(3) * cos(t17) + t20;
t32 = t34 * t13;
t31 = cos(pkin(10)) * pkin(2) + pkin(1) + t33;
t28 = r_i_i_C(1) + pkin(8) + pkin(7) + qJ(2);
t18 = sin(qJ(1));
t27 = t18 * t12;
t19 = cos(qJ(1));
t26 = t19 * t12;
t23 = r_i_i_C(2) * t27 + t32 * t18;
t22 = r_i_i_C(2) * t26 + t32 * t19;
t21 = -pkin(3) * sin(t17) - pkin(4) * t12;
t1 = [-t31 * t18 + t28 * t19, t18, t21 * t19 + t22, -pkin(4) * t26 + t22, t26, 0; t28 * t18 + t31 * t19, -t19, t21 * t18 + t23, -pkin(4) * t27 + t23, t27, 0; 0, 0, t33, t20, -t13, 0;];
Ja_transl  = t1;
