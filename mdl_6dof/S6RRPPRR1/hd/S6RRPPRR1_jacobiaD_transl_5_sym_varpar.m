% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:28:06
% EndTime: 2019-02-26 21:28:07
% DurationCPUTime: 0.19s
% Computational Cost: add. (260->40), mult. (402->58), div. (0->0), fcn. (344->8), ass. (0->34)
t70 = qJ(2) + pkin(10);
t68 = sin(t70);
t105 = pkin(3) + pkin(4);
t69 = cos(t70);
t81 = t105 * t68 + sin(qJ(2)) * pkin(2) - qJ(4) * t69;
t78 = -t81 * qJD(2) + t68 * qJD(4);
t111 = t78 - (pkin(8) + r_i_i_C(3) - qJ(3) - pkin(7)) * qJD(1);
t72 = sin(qJ(5));
t75 = cos(qJ(5));
t110 = -t68 * t75 + t69 * t72;
t107 = qJD(2) - qJD(5);
t86 = t68 * t72 + t69 * t75;
t62 = t107 * t86;
t109 = -qJ(4) * t68 - t105 * t69 - cos(qJ(2)) * pkin(2);
t99 = qJD(2) * t68;
t106 = t110 * qJD(5) + t75 * t99;
t74 = sin(qJ(1));
t101 = qJD(1) * t74;
t77 = cos(qJ(1));
t100 = qJD(1) * t77;
t98 = qJD(2) * t69;
t92 = t77 * t98;
t84 = qJD(1) * t110;
t57 = t62 * t77 + t74 * t84;
t83 = qJD(1) * t86;
t58 = t106 * t77 - t72 * t92 + t74 * t83;
t91 = t57 * r_i_i_C(1) + t58 * r_i_i_C(2);
t59 = -t74 * t62 + t77 * t84;
t60 = t77 * t83 + (t72 * t98 - t106) * t74;
t90 = -t59 * r_i_i_C(1) - t60 * r_i_i_C(2);
t89 = -r_i_i_C(1) * t107 * t110 - t62 * r_i_i_C(2);
t80 = qJD(3) + (-pkin(1) + t109) * qJD(1);
t79 = t109 * qJD(2) + qJD(4) * t69;
t1 = [-t60 * r_i_i_C(1) + t59 * r_i_i_C(2) - t111 * t74 + t80 * t77, t81 * t101 + t79 * t77 - t91, t100, -t68 * t101 + t92, t91, 0; -t58 * r_i_i_C(1) + t57 * r_i_i_C(2) + t111 * t77 + t80 * t74, -t100 * t81 + t79 * t74 - t90, t101, t68 * t100 + t74 * t98, t90, 0; 0, t78 - t89, 0, t99, t89, 0;];
JaD_transl  = t1;
