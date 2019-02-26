% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR13
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR13_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR13_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobia_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:55:38
% EndTime: 2019-02-26 20:55:38
% DurationCPUTime: 0.12s
% Computational Cost: add. (118->34), mult. (317->57), div. (0->0), fcn. (413->10), ass. (0->34)
t44 = r_i_i_C(1) + pkin(9);
t21 = cos(pkin(6));
t16 = sin(pkin(12));
t25 = cos(qJ(1));
t34 = t25 * t16;
t19 = cos(pkin(12));
t23 = sin(qJ(1));
t35 = t23 * t19;
t11 = t21 * t34 + t35;
t22 = sin(qJ(3));
t24 = cos(qJ(3));
t32 = t25 * t19;
t37 = t23 * t16;
t10 = -t21 * t32 + t37;
t17 = sin(pkin(7));
t20 = cos(pkin(7));
t18 = sin(pkin(6));
t33 = t25 * t18;
t27 = t10 * t20 + t17 * t33;
t43 = -t11 * t24 + t27 * t22;
t1 = t11 * t22 + t27 * t24;
t42 = -r_i_i_C(2) + pkin(3);
t39 = t17 * t21;
t38 = t20 * t24;
t36 = t23 * t18;
t31 = r_i_i_C(3) + qJ(4);
t30 = t18 * qJ(2);
t29 = t17 * t36;
t13 = -t21 * t37 + t32;
t12 = -t21 * t35 - t34;
t7 = -t24 * t39 + (t16 * t22 - t19 * t38) * t18;
t6 = t13 * t24 + (t12 * t20 + t29) * t22;
t5 = -t12 * t38 + t13 * t22 - t24 * t29;
t2 = [-t11 * pkin(2) - t23 * pkin(1) + t25 * t30 + t42 * t43 - t31 * t1 + t44 * (-t10 * t17 + t20 * t33) t36, t31 * t6 - t42 * t5, t5, 0, 0; t25 * pkin(1) + t13 * pkin(2) + t23 * t30 + t31 * t5 + t42 * t6 + t44 * (-t12 * t17 + t20 * t36) -t33, -t42 * t1 - t31 * t43, t1, 0, 0; 0, t21, t31 * (t22 * t39 + (t19 * t20 * t22 + t16 * t24) * t18) - t42 * t7, t7, 0, 0;];
Ja_transl  = t2;
