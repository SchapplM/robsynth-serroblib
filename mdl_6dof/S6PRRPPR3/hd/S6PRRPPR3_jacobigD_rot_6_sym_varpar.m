% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPPR3_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_jacobigD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:59:24
% EndTime: 2019-02-26 19:59:24
% DurationCPUTime: 0.06s
% Computational Cost: add. (12->10), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
t127 = sin(pkin(6));
t131 = sin(qJ(2));
t140 = t127 * t131;
t132 = cos(qJ(3));
t139 = t127 * t132;
t129 = cos(pkin(6));
t138 = t129 * t131;
t133 = cos(qJ(2));
t137 = t129 * t133;
t136 = qJD(2) * t132;
t126 = sin(pkin(10));
t128 = cos(pkin(10));
t135 = t126 * t133 + t128 * t138;
t134 = -t126 * t138 + t128 * t133;
t130 = sin(qJ(3));
t1 = [0, 0, t134 * qJD(2), 0, 0 (t126 * t139 - t134 * t130) * qJD(3) + (-t126 * t137 - t128 * t131) * t136; 0, 0, t135 * qJD(2), 0, 0 (-t128 * t139 - t135 * t130) * qJD(3) + (-t126 * t131 + t128 * t137) * t136; 0, 0, qJD(2) * t140, 0, 0, t127 * t133 * t136 + (t129 * t132 - t130 * t140) * qJD(3);];
JgD_rot  = t1;
