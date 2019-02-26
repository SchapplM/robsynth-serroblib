% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function Ja_transl = S6RRRRPP3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:26:19
% EndTime: 2019-02-26 22:26:20
% DurationCPUTime: 0.16s
% Computational Cost: add. (204->29), mult. (251->40), div. (0->0), fcn. (281->8), ass. (0->30)
t27 = sin(qJ(4));
t30 = cos(qJ(4));
t38 = pkin(4) + r_i_i_C(3) + qJ(6);
t41 = r_i_i_C(2) + qJ(5);
t56 = t41 * t27 + t38 * t30;
t53 = pkin(5) + pkin(9) + r_i_i_C(1);
t26 = qJ(2) + qJ(3);
t24 = cos(t26);
t52 = t24 * t53;
t23 = sin(t26);
t51 = t24 * pkin(3) + t53 * t23;
t25 = cos(qJ(2)) * pkin(2);
t50 = pkin(1) + t25 + t51;
t29 = sin(qJ(1));
t45 = t29 * t27;
t44 = t29 * t30;
t31 = cos(qJ(1));
t43 = t31 * t27;
t42 = t31 * t30;
t39 = t29 * t52;
t37 = t31 * t52;
t35 = t56 * t24 + t51;
t34 = (-pkin(3) - t56) * t23;
t33 = -sin(qJ(2)) * pkin(2) + t34;
t32 = -pkin(8) - pkin(7);
t4 = t24 * t42 + t45;
t3 = t24 * t43 - t44;
t2 = t24 * t44 - t43;
t1 = t24 * t45 + t42;
t5 = [-t41 * t1 - t38 * t2 - t50 * t29 - t31 * t32, t33 * t31 + t37, t31 * t34 + t37, -t38 * t3 + t41 * t4, t3, t4; -t29 * t32 + t41 * t3 + t50 * t31 + t38 * t4, t33 * t29 + t39, t29 * t34 + t39, -t38 * t1 + t41 * t2, t1, t2; 0, t25 + t35, t35 (-t38 * t27 + t41 * t30) * t23, t23 * t27, t23 * t30;];
Ja_transl  = t5;
