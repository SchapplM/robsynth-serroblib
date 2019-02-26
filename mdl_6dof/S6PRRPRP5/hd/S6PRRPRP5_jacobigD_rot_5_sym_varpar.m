% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPRP5_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_jacobigD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:03:31
% EndTime: 2019-02-26 20:03:31
% DurationCPUTime: 0.07s
% Computational Cost: add. (12->10), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
t128 = sin(pkin(6));
t132 = sin(qJ(2));
t141 = t128 * t132;
t133 = cos(qJ(3));
t140 = t128 * t133;
t130 = cos(pkin(6));
t139 = t130 * t132;
t134 = cos(qJ(2));
t138 = t130 * t134;
t137 = qJD(2) * t133;
t127 = sin(pkin(10));
t129 = cos(pkin(10));
t136 = t127 * t134 + t129 * t139;
t135 = -t127 * t139 + t129 * t134;
t131 = sin(qJ(3));
t1 = [0, 0, t135 * qJD(2), 0 (t127 * t140 - t135 * t131) * qJD(3) + (-t127 * t138 - t129 * t132) * t137, 0; 0, 0, t136 * qJD(2), 0 (-t129 * t140 - t136 * t131) * qJD(3) + (-t127 * t132 + t129 * t138) * t137, 0; 0, 0, qJD(2) * t141, 0, t128 * t134 * t137 + (t130 * t133 - t131 * t141) * qJD(3), 0;];
JgD_rot  = t1;
