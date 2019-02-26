% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPP5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:27:38
% EndTime: 2019-02-26 22:27:38
% DurationCPUTime: 0.14s
% Computational Cost: add. (209->33), mult. (250->46), div. (0->0), fcn. (288->8), ass. (0->31)
t38 = r_i_i_C(2) + qJ(5);
t24 = cos(qJ(3));
t17 = t24 * pkin(3) + pkin(2);
t25 = cos(qJ(2));
t22 = sin(qJ(2));
t35 = r_i_i_C(3) + qJ(6) - pkin(9) - pkin(8);
t31 = t35 * t22;
t45 = -t25 * t17 - pkin(1) + t31;
t36 = pkin(4) + pkin(5) + r_i_i_C(1);
t20 = qJ(3) + qJ(4);
t18 = sin(t20);
t19 = cos(t20);
t44 = t38 * t18 + t36 * t19 + t17;
t21 = sin(qJ(3));
t43 = pkin(3) * t21;
t23 = sin(qJ(1));
t41 = t23 * t25;
t26 = cos(qJ(1));
t40 = t26 * t18;
t39 = t26 * t19;
t37 = t38 * t19 * t22;
t34 = pkin(7) + t43;
t32 = t36 * t18;
t10 = t19 * t41 - t40;
t9 = t18 * t41 + t39;
t30 = t38 * t10 - t36 * t9;
t11 = -t23 * t19 + t25 * t40;
t12 = t23 * t18 + t25 * t39;
t29 = -t36 * t11 + t38 * t12;
t28 = -t44 * t22 - t35 * t25;
t1 = [-t36 * t10 + t45 * t23 + t34 * t26 - t38 * t9, t28 * t26 (-t21 * t25 * t26 + t23 * t24) * pkin(3) + t29, t29, t11, -t26 * t22; t38 * t11 + t36 * t12 + t34 * t23 - t45 * t26, t28 * t23 (-t21 * t41 - t24 * t26) * pkin(3) + t30, t30, t9, -t23 * t22; 0, t44 * t25 - t31 (-t32 - t43) * t22 + t37, -t22 * t32 + t37, t22 * t18, t25;];
Ja_transl  = t1;
