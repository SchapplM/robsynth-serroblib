% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRP3_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP3_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_jacobia_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:10:35
% EndTime: 2019-02-26 22:10:35
% DurationCPUTime: 0.10s
% Computational Cost: add. (91->19), mult. (103->23), div. (0->0), fcn. (110->8), ass. (0->20)
t40 = r_i_i_C(3) + qJ(4);
t15 = qJ(2) + qJ(3);
t12 = sin(t15);
t13 = cos(t15);
t17 = cos(pkin(10));
t27 = -r_i_i_C(1) * t17 - pkin(3);
t16 = sin(pkin(10));
t34 = r_i_i_C(2) * t16;
t23 = t40 * t12 + (-t27 - t34) * t13;
t39 = cos(qJ(2)) * pkin(2) + t23;
t38 = t12 * t34 + t40 * t13;
t36 = pkin(1) + t39;
t19 = sin(qJ(1));
t30 = t38 * t19;
t20 = cos(qJ(1));
t29 = t38 * t20;
t26 = t27 * t12;
t24 = r_i_i_C(1) * t16 + r_i_i_C(2) * t17 + pkin(7) + pkin(8);
t22 = -sin(qJ(2)) * pkin(2) + t26;
t1 = [-t36 * t19 + t24 * t20, t22 * t20 + t29, t20 * t26 + t29, t20 * t12, 0, 0; t24 * t19 + t36 * t20, t22 * t19 + t30, t19 * t26 + t30, t19 * t12, 0, 0; 0, t39, t23, -t13, 0, 0;];
Ja_transl  = t1;
