% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_transl = S6RRRRPP5_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP5_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_jacobia_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:27:43
% EndTime: 2019-02-26 22:27:43
% DurationCPUTime: 0.12s
% Computational Cost: add. (168->31), mult. (203->44), div. (0->0), fcn. (232->8), ass. (0->31)
t41 = pkin(4) + r_i_i_C(1);
t34 = r_i_i_C(3) + qJ(5);
t22 = cos(qJ(3));
t15 = t22 * pkin(3) + pkin(2);
t23 = cos(qJ(2));
t20 = sin(qJ(2));
t39 = r_i_i_C(2) + pkin(9) + pkin(8);
t30 = t39 * t20;
t43 = t23 * t15 + pkin(1) + t30;
t18 = qJ(3) + qJ(4);
t16 = sin(t18);
t17 = cos(t18);
t42 = t34 * t16 + t41 * t17 + t15;
t19 = sin(qJ(3));
t40 = pkin(3) * t19;
t21 = sin(qJ(1));
t37 = t21 * t23;
t24 = cos(qJ(1));
t36 = t24 * t16;
t35 = t24 * t17;
t33 = t34 * t17 * t20;
t32 = t41 * t16;
t31 = pkin(7) + t40;
t7 = t16 * t37 + t35;
t8 = t17 * t37 - t36;
t28 = t34 * t8 - t41 * t7;
t10 = t21 * t16 + t23 * t35;
t9 = -t21 * t17 + t23 * t36;
t27 = t34 * t10 - t41 * t9;
t26 = -t42 * t20 + t39 * t23;
t1 = [-t43 * t21 + t31 * t24 - t34 * t7 - t41 * t8, t26 * t24 (-t19 * t23 * t24 + t21 * t22) * pkin(3) + t27, t27, t9, 0; t41 * t10 + t31 * t21 + t43 * t24 + t34 * t9, t26 * t21 (-t19 * t37 - t22 * t24) * pkin(3) + t28, t28, t7, 0; 0, t42 * t23 + t30 (-t32 - t40) * t20 + t33, -t20 * t32 + t33, t20 * t16, 0;];
Ja_transl  = t1;
