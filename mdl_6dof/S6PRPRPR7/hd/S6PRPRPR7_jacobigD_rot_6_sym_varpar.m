% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRPR7_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_jacobigD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:49:47
% EndTime: 2019-02-26 19:49:47
% DurationCPUTime: 0.04s
% Computational Cost: add. (12->11), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
t126 = sin(pkin(6));
t129 = sin(qJ(4));
t139 = t126 * t129;
t132 = cos(qJ(2));
t138 = t126 * t132;
t128 = cos(pkin(6));
t130 = sin(qJ(2));
t137 = t128 * t130;
t136 = t128 * t132;
t135 = qJD(2) * t129;
t125 = sin(pkin(10));
t127 = cos(pkin(10));
t134 = -t125 * t130 + t127 * t136;
t133 = t125 * t136 + t127 * t130;
t131 = cos(qJ(4));
t1 = [0, 0, 0, -t133 * qJD(2), 0 (-t125 * t139 + t133 * t131) * qJD(4) + (-t125 * t137 + t127 * t132) * t135; 0, 0, 0, t134 * qJD(2), 0 (t127 * t139 - t134 * t131) * qJD(4) + (t125 * t132 + t127 * t137) * t135; 0, 0, 0, qJD(2) * t138, 0, t126 * t130 * t135 + (-t128 * t129 - t131 * t138) * qJD(4);];
JgD_rot  = t1;
