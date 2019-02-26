% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP7_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP7_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_jacobiaD_transl_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:49:25
% EndTime: 2019-02-26 21:49:25
% DurationCPUTime: 0.18s
% Computational Cost: add. (122->38), mult. (376->64), div. (0->0), fcn. (325->6), ass. (0->33)
t59 = sin(qJ(4));
t60 = sin(qJ(2));
t62 = cos(qJ(4));
t63 = cos(qJ(2));
t72 = t59 * t63 - t60 * t62;
t90 = pkin(2) + pkin(3);
t91 = -qJ(3) * t63 + t90 * t60;
t94 = t91 * qJD(2) - t60 * qJD(3);
t71 = t59 * t60 + t62 * t63;
t93 = qJD(2) - qJD(4);
t54 = t93 * t71;
t84 = qJD(2) * t60;
t92 = t72 * qJD(4) + t62 * t84;
t61 = sin(qJ(1));
t86 = qJD(1) * t61;
t64 = cos(qJ(1));
t85 = qJD(1) * t64;
t83 = qJD(2) * t63;
t82 = qJD(2) * t64;
t80 = pkin(7) - pkin(8) - r_i_i_C(3);
t76 = t63 * t82;
t68 = qJD(1) * t72;
t49 = t54 * t64 + t61 * t68;
t67 = qJD(1) * t71;
t50 = -t59 * t76 + t61 * t67 + t92 * t64;
t75 = t49 * r_i_i_C(1) + t50 * r_i_i_C(2);
t51 = -t61 * t54 + t64 * t68;
t52 = t64 * t67 + (t59 * t83 - t92) * t61;
t74 = -t51 * r_i_i_C(1) - t52 * r_i_i_C(2);
t73 = -t93 * t72 * r_i_i_C(1) - t54 * r_i_i_C(2);
t70 = -qJ(3) * t60 - t90 * t63;
t66 = -pkin(1) + t70;
t1 = [-t52 * r_i_i_C(1) + t51 * r_i_i_C(2) + t94 * t61 + (-t80 * t61 + t66 * t64) * qJD(1) (-qJ(3) * t82 + t90 * t86) * t60 + (-qJ(3) * t86 + (-t90 * qJD(2) + qJD(3)) * t64) * t63 - t75, -t60 * t86 + t76, t75, 0, 0; -t50 * r_i_i_C(1) + t49 * r_i_i_C(2) - t94 * t64 + (t66 * t61 + t80 * t64) * qJD(1), -t91 * t85 + (qJD(2) * t70 + qJD(3) * t63) * t61 - t74, t60 * t85 + t61 * t83, t74, 0, 0; 0, -t94 - t73, t84, t73, 0, 0;];
JaD_transl  = t1;
