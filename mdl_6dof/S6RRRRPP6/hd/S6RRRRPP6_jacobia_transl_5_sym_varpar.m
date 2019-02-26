% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP6
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
% Datum: 2019-02-26 22:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPP6_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP6_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_jacobia_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:28:23
% EndTime: 2019-02-26 22:28:23
% DurationCPUTime: 0.12s
% Computational Cost: add. (168->32), mult. (203->46), div. (0->0), fcn. (232->8), ass. (0->31)
t42 = pkin(4) - r_i_i_C(2);
t34 = r_i_i_C(3) + qJ(5);
t23 = cos(qJ(3));
t16 = t23 * pkin(3) + pkin(2);
t24 = cos(qJ(2));
t21 = sin(qJ(2));
t40 = r_i_i_C(1) + pkin(9) + pkin(8);
t31 = t40 * t21;
t44 = t24 * t16 + pkin(1) + t31;
t19 = qJ(3) + qJ(4);
t17 = sin(t19);
t18 = cos(t19);
t43 = t34 * t17 + t42 * t18 + t16;
t20 = sin(qJ(3));
t41 = pkin(3) * t20;
t38 = t21 * t17;
t22 = sin(qJ(1));
t37 = t22 * t24;
t25 = cos(qJ(1));
t36 = t25 * t17;
t35 = t25 * t18;
t33 = t34 * t18 * t21 + r_i_i_C(2) * t38;
t32 = pkin(7) + t41;
t7 = t17 * t37 + t35;
t8 = t18 * t37 - t36;
t29 = t34 * t8 - t42 * t7;
t10 = t22 * t17 + t24 * t35;
t9 = -t22 * t18 + t24 * t36;
t28 = t34 * t10 - t42 * t9;
t27 = -t43 * t21 + t40 * t24;
t1 = [-t44 * t22 + t32 * t25 - t34 * t7 - t42 * t8, t27 * t25 (-t20 * t24 * t25 + t22 * t23) * pkin(3) + t28, t28, t9, 0; t42 * t10 + t32 * t22 + t44 * t25 + t34 * t9, t27 * t22 (-t20 * t37 - t23 * t25) * pkin(3) + t29, t29, t7, 0; 0, t43 * t24 + t31 (-pkin(4) * t17 - t41) * t21 + t33, -pkin(4) * t38 + t33, t38, 0;];
Ja_transl  = t1;
