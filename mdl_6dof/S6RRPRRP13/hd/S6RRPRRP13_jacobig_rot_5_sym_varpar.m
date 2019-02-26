% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRP13_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_jacobig_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:52:49
% EndTime: 2019-02-26 21:52:49
% DurationCPUTime: 0.03s
% Computational Cost: add. (8->8), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->15)
t100 = sin(pkin(6));
t104 = sin(qJ(1));
t113 = t104 * t100;
t103 = sin(qJ(2));
t112 = t104 * t103;
t106 = cos(qJ(2));
t111 = t104 * t106;
t107 = cos(qJ(1));
t110 = t107 * t100;
t109 = t107 * t103;
t108 = t107 * t106;
t105 = cos(qJ(4));
t102 = sin(qJ(4));
t101 = cos(pkin(6));
t1 = [0, t113, 0, -t101 * t112 + t108, t102 * t113 - (t101 * t111 + t109) * t105, 0; 0, -t110, 0, t101 * t109 + t111, -t102 * t110 - (-t101 * t108 + t112) * t105, 0; 1, t101, 0, t100 * t103, t100 * t106 * t105 + t101 * t102, 0;];
Jg_rot  = t1;
