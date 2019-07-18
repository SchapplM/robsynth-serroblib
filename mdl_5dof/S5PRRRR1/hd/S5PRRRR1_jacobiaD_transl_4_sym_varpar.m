% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRRR1_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_jacobiaD_transl_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_jacobiaD_transl_4_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR1_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_jacobiaD_transl_4_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:29:14
% EndTime: 2019-07-18 13:29:14
% DurationCPUTime: 0.10s
% Computational Cost: add. (75->24), mult. (110->40), div. (0->0), fcn. (71->6), ass. (0->27)
t93 = sin(qJ(3));
t112 = pkin(2) * t93;
t100 = qJD(3) * t112;
t91 = qJD(3) + qJD(4);
t92 = qJ(3) + qJ(4);
t90 = cos(t92);
t110 = r_i_i_C(2) * t90;
t89 = sin(t92);
t111 = r_i_i_C(1) * t89;
t99 = t110 + t111;
t113 = t99 * t91 + t100;
t109 = t89 * t91;
t108 = t90 * t91;
t107 = r_i_i_C(1) * t109 + r_i_i_C(2) * t108;
t94 = sin(qJ(2));
t106 = qJD(2) * t94;
t96 = cos(qJ(2));
t105 = qJD(2) * t96;
t95 = cos(qJ(3));
t104 = qJD(3) * t95;
t103 = r_i_i_C(1) * t108;
t102 = r_i_i_C(2) * t109;
t101 = qJD(2) * t110;
t98 = -pkin(2) * t95 - r_i_i_C(1) * t90 + r_i_i_C(2) * t89;
t97 = t94 * t101 + t106 * t111 + (t102 - t103) * t96;
t83 = t94 * t102;
t1 = [0, t113 * t94 + (-r_i_i_C(3) * t94 + t98 * t96) * qJD(2), (-t96 * t104 + t93 * t106) * pkin(2) + t97, t97, 0; 0, 0, t100 + t107, t107, 0; 0, -t113 * t96 + (r_i_i_C(3) * t96 + t98 * t94) * qJD(2), t83 + (-pkin(2) * t104 - t103) * t94 + (-t99 - t112) * t105, -t96 * t101 + t83 + (-t89 * t105 - t94 * t108) * r_i_i_C(1), 0;];
JaD_transl  = t1;
