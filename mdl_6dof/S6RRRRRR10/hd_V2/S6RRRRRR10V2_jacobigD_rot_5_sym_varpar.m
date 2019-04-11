% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRR10V2_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobigD_rot_5_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:56:32
% EndTime: 2019-04-11 14:56:32
% DurationCPUTime: 0.03s
% Computational Cost: add. (30->11), mult. (43->24), div. (0->0), fcn. (43->6), ass. (0->15)
t157 = qJD(2) + qJD(3);
t160 = sin(qJ(1));
t166 = t157 * t160;
t162 = cos(qJ(1));
t165 = t157 * t162;
t153 = qJD(1) * t160;
t154 = qJD(1) * t162;
t158 = qJ(2) + qJ(3);
t156 = cos(t158);
t164 = qJD(1) * t156 - qJD(4);
t161 = cos(qJ(4));
t163 = (qJD(4) * t156 - qJD(1)) * t161;
t159 = sin(qJ(4));
t155 = sin(t158);
t1 = [0, t154, t154, -t155 * t153 + t156 * t165, t162 * t163 + (-t155 * t165 - t164 * t160) * t159, 0; 0, t153, t153, t155 * t154 + t156 * t166, t160 * t163 + (-t155 * t166 + t164 * t162) * t159, 0; 0, 0, 0, t157 * t155, t155 * qJD(4) * t161 + t157 * t156 * t159, 0;];
JgD_rot  = t1;
