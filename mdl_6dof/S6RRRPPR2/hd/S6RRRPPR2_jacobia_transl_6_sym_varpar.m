% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:03:55
% EndTime: 2019-02-26 22:03:56
% DurationCPUTime: 0.13s
% Computational Cost: add. (198->32), mult. (149->42), div. (0->0), fcn. (162->10), ass. (0->28)
t27 = sin(qJ(6));
t29 = cos(qJ(6));
t54 = r_i_i_C(1) * t27 + r_i_i_C(2) * t29;
t26 = qJ(2) + qJ(3);
t22 = pkin(10) + t26;
t19 = sin(t22);
t20 = cos(t22);
t40 = pkin(4) + pkin(9) + r_i_i_C(3);
t53 = t40 * t20 + t19 * qJ(5) + pkin(3) * cos(t26);
t51 = (qJ(5) + t54) * t20;
t31 = -t40 * t19 - pkin(3) * sin(t26);
t24 = cos(qJ(2)) * pkin(2);
t49 = pkin(1) + t24 + t53;
t45 = pkin(5) + qJ(4) + pkin(8) + pkin(7);
t28 = sin(qJ(1));
t44 = t28 * t27;
t43 = t28 * t29;
t30 = cos(qJ(1));
t42 = t30 * t19;
t39 = t51 * t28;
t38 = t51 * t30;
t33 = -sin(qJ(2)) * pkin(2) + t31;
t32 = t54 * t19 + t53;
t4 = -t19 * t44 + t29 * t30;
t3 = t19 * t43 + t27 * t30;
t2 = t27 * t42 + t43;
t1 = t29 * t42 - t44;
t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t49 * t28 + t45 * t30, t33 * t30 + t38, t31 * t30 + t38, t28, t42, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t45 * t28 + t49 * t30, t33 * t28 + t39, t31 * t28 + t39, -t30, t28 * t19, r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, t24 + t32, t32, 0, -t20 (-r_i_i_C(1) * t29 + r_i_i_C(2) * t27) * t20;];
Ja_transl  = t5;
