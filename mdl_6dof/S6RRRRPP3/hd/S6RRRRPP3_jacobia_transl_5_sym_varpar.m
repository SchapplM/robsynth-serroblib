% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP3
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
% Datum: 2019-02-26 22:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPP3_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP3_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_jacobia_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:26:19
% EndTime: 2019-02-26 22:26:19
% DurationCPUTime: 0.13s
% Computational Cost: add. (157->31), mult. (196->42), div. (0->0), fcn. (217->8), ass. (0->30)
t23 = qJ(2) + qJ(3);
t20 = sin(t23);
t21 = cos(t23);
t27 = cos(qJ(4));
t52 = pkin(9) + r_i_i_C(1);
t53 = r_i_i_C(2) * t20 * t27 + t21 * t52;
t47 = pkin(4) - r_i_i_C(2);
t24 = sin(qJ(4));
t37 = r_i_i_C(3) + qJ(5);
t51 = t37 * t24;
t49 = t21 * pkin(3) + t52 * t20;
t22 = cos(qJ(2)) * pkin(2);
t48 = pkin(1) + t22 + t49;
t26 = sin(qJ(1));
t41 = t24 * t26;
t40 = t26 * t27;
t28 = cos(qJ(1));
t39 = t27 * t28;
t38 = t28 * t24;
t36 = t53 * t26;
t34 = t53 * t28;
t32 = (-pkin(4) * t27 - pkin(3) - t51) * t20;
t31 = t49 + (t27 * t47 + t51) * t21;
t30 = -sin(qJ(2)) * pkin(2) + t32;
t29 = -pkin(8) - pkin(7);
t4 = t21 * t39 + t41;
t3 = t21 * t38 - t40;
t2 = t21 * t40 - t38;
t1 = t21 * t41 + t39;
t5 = [-t37 * t1 - t47 * t2 - t48 * t26 - t28 * t29, t30 * t28 + t34, t28 * t32 + t34, -t47 * t3 + t37 * t4, t3, 0; -t26 * t29 + t48 * t28 + t37 * t3 + t47 * t4, t30 * t26 + t36, t26 * t32 + t36, -t47 * t1 + t37 * t2, t1, 0; 0, t22 + t31, t31 (-t47 * t24 + t37 * t27) * t20, t20 * t24, 0;];
Ja_transl  = t5;
