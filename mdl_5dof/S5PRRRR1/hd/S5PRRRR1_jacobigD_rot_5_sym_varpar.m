% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
%
% Output:
% JgD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S5PRRRR1_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_jacobigD_rot_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_jacobigD_rot_5_sym_varpar: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_jacobigD_rot_5_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:29:14
% EndTime: 2019-07-18 13:29:14
% DurationCPUTime: 0.02s
% Computational Cost: add. (11->5), mult. (13->8), div. (0->0), fcn. (13->4), ass. (0->9)
t109 = qJD(3) + qJD(4);
t110 = qJ(3) + qJ(4);
t113 = cos(t110) * t109;
t111 = sin(qJ(2));
t105 = qJD(2) * t111;
t112 = cos(qJ(2));
t106 = qJD(2) * t112;
t107 = sin(t110);
t1 = [0, 0, t106, t106, -t107 * t105 + t112 * t113; 0, 0, 0, 0, -t109 * t107; 0, 0, t105, t105, t107 * t106 + t111 * t113;];
JgD_rot  = t1;
