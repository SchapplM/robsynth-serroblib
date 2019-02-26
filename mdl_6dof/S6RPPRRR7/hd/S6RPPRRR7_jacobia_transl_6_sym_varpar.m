% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRR7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:38:11
% EndTime: 2019-02-26 20:38:11
% DurationCPUTime: 0.13s
% Computational Cost: add. (158->35), mult. (127->43), div. (0->0), fcn. (139->9), ass. (0->31)
t20 = pkin(10) + qJ(4);
t18 = qJ(5) + t20;
t14 = sin(t18);
t41 = pkin(9) + r_i_i_C(3);
t43 = t41 * t14;
t15 = cos(t18);
t42 = t41 * t15;
t39 = pkin(4) * sin(t20);
t38 = pkin(4) * cos(t20);
t21 = sin(qJ(6));
t37 = r_i_i_C(2) * t21;
t36 = pkin(1) + pkin(8) + pkin(7) + qJ(3);
t24 = cos(qJ(1));
t34 = t21 * t24;
t22 = sin(qJ(1));
t33 = t22 * t21;
t23 = cos(qJ(6));
t32 = t22 * t23;
t31 = t23 * t24;
t30 = t15 * t37;
t29 = -r_i_i_C(1) * t23 - pkin(5);
t28 = t22 * t43 + (pkin(5) * t22 + r_i_i_C(1) * t32) * t15;
t27 = t42 + (t29 + t37) * t14;
t26 = pkin(5) * t14 - t42 + qJ(2) + t39 + sin(pkin(10)) * pkin(3);
t25 = t29 * t15 - t43;
t6 = t24 * t30;
t4 = t14 * t31 - t33;
t3 = t14 * t34 + t32;
t2 = t14 * t32 + t34;
t1 = -t14 * t33 + t31;
t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t36 * t22 + t26 * t24, t22, t24 (-t30 + t38) * t22 + t28, -t22 * t30 + t28, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t26 * t22 + t36 * t24, -t24, t22, t6 + (t25 - t38) * t24, t25 * t24 + t6, r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, 0, 0, t27 - t39, t27 (-r_i_i_C(1) * t21 - r_i_i_C(2) * t23) * t15;];
Ja_transl  = t5;
