% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRRR14V3
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR14V3_jacobiaD_transl_3_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobiaD_transl_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_jacobiaD_transl_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14V3_jacobiaD_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobiaD_transl_3_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:12:05
% EndTime: 2019-04-12 15:12:05
% DurationCPUTime: 0.07s
% Computational Cost: add. (31->15), mult. (100->32), div. (0->0), fcn. (71->4), ass. (0->14)
t132 = sin(qJ(2));
t134 = cos(qJ(2));
t145 = r_i_i_C(3) + qJ(3);
t139 = -r_i_i_C(1) * t132 + t145 * t134;
t136 = t139 * qJD(2) + t132 * qJD(3);
t148 = r_i_i_C(2) * qJD(1) + t136;
t133 = sin(qJ(1));
t143 = qJD(1) * t133;
t135 = cos(qJ(1));
t142 = qJD(1) * t135;
t141 = qJD(2) * t135;
t138 = -r_i_i_C(1) * t134 - t145 * t132;
t137 = qJD(1) * t138;
t1 = [-t148 * t133 + t135 * t137 (r_i_i_C(1) * t143 - t145 * t141) * t132 + ((-r_i_i_C(1) * qJD(2) + qJD(3)) * t135 - t145 * t143) * t134, -t132 * t143 + t134 * t141, 0, 0, 0; t133 * t137 + t148 * t135, t139 * t142 + (t138 * qJD(2) + qJD(3) * t134) * t133, t133 * qJD(2) * t134 + t132 * t142, 0, 0, 0; 0, t136, qJD(2) * t132, 0, 0, 0;];
JaD_transl  = t1;
