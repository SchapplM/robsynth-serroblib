% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function Ja_transl = S6RRRRPP6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:28:13
% EndTime: 2019-02-26 22:28:14
% DurationCPUTime: 0.13s
% Computational Cost: add. (214->31), mult. (255->45), div. (0->0), fcn. (295->8), ass. (0->32)
t36 = r_i_i_C(2) + qJ(5);
t22 = cos(qJ(3));
t15 = t22 * pkin(3) + pkin(2);
t23 = cos(qJ(2));
t20 = sin(qJ(2));
t34 = pkin(5) + r_i_i_C(1) + pkin(9) + pkin(8);
t30 = t34 * t20;
t43 = t23 * t15 + pkin(1) + t30;
t33 = pkin(4) + r_i_i_C(3) + qJ(6);
t18 = qJ(3) + qJ(4);
t16 = sin(t18);
t17 = cos(t18);
t42 = t36 * t16 + t33 * t17 + t15;
t19 = sin(qJ(3));
t41 = pkin(3) * t19;
t40 = t20 * t17;
t21 = sin(qJ(1));
t39 = t21 * t23;
t24 = cos(qJ(1));
t38 = t24 * t16;
t37 = t24 * t17;
t35 = t36 * t40;
t32 = pkin(7) + t41;
t29 = t33 * t16;
t7 = t16 * t39 + t37;
t8 = t17 * t39 - t38;
t28 = -t33 * t7 + t36 * t8;
t10 = t21 * t16 + t23 * t37;
t9 = -t21 * t17 + t23 * t38;
t27 = t36 * t10 - t33 * t9;
t26 = -t42 * t20 + t34 * t23;
t1 = [-t43 * t21 + t32 * t24 - t33 * t8 - t36 * t7, t26 * t24 (-t19 * t23 * t24 + t21 * t22) * pkin(3) + t27, t27, t9, t10; t33 * t10 + t32 * t21 + t43 * t24 + t36 * t9, t26 * t21 (-t19 * t39 - t22 * t24) * pkin(3) + t28, t28, t7, t8; 0, t42 * t23 + t30 (-t29 - t41) * t20 + t35, -t20 * t29 + t35, t20 * t16, t40;];
Ja_transl  = t1;
