% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:52:12
% EndTime: 2019-02-26 20:52:12
% DurationCPUTime: 0.13s
% Computational Cost: add. (161->36), mult. (130->44), div. (0->0), fcn. (142->10), ass. (0->31)
t20 = qJ(3) + pkin(10);
t18 = qJ(5) + t20;
t15 = sin(t18);
t40 = pkin(9) + r_i_i_C(3);
t42 = t40 * t15;
t16 = cos(t18);
t41 = t40 * t16;
t21 = sin(qJ(6));
t38 = r_i_i_C(2) * t21;
t37 = pkin(1) + pkin(8) + qJ(4) + pkin(7);
t25 = cos(qJ(1));
t35 = t21 * t25;
t23 = sin(qJ(1));
t34 = t23 * t21;
t24 = cos(qJ(6));
t33 = t23 * t24;
t32 = t24 * t25;
t31 = t16 * t38;
t30 = -r_i_i_C(1) * t24 - pkin(5);
t29 = t23 * t42 + (pkin(5) * t23 + r_i_i_C(1) * t33) * t16;
t7 = sin(qJ(3)) * pkin(3) + pkin(4) * sin(t20);
t28 = t41 + (t30 + t38) * t15;
t27 = pkin(5) * t15 + qJ(2) - t41 + t7;
t26 = t30 * t16 - t42;
t8 = pkin(4) * cos(t20) + cos(qJ(3)) * pkin(3);
t6 = t25 * t31;
t4 = t15 * t32 - t34;
t3 = t15 * t35 + t33;
t2 = t15 * t33 + t35;
t1 = -t15 * t34 + t32;
t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t37 * t23 + t27 * t25, t23 (t8 - t31) * t23 + t29, t25, -t23 * t31 + t29, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t27 * t23 + t37 * t25, -t25, t6 + (-t8 + t26) * t25, t23, t26 * t25 + t6, r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, 0, t28 - t7, 0, t28 (-r_i_i_C(1) * t21 - r_i_i_C(2) * t24) * t16;];
Ja_transl  = t5;
