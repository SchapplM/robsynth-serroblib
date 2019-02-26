% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRP3
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPRP3_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_jacobig_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:02:32
% EndTime: 2019-02-26 20:02:32
% DurationCPUTime: 0.02s
% Computational Cost: add. (9->9), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->13)
t90 = sin(pkin(10));
t91 = sin(pkin(6));
t101 = t90 * t91;
t92 = cos(pkin(10));
t100 = t92 * t91;
t93 = cos(pkin(6));
t95 = sin(qJ(2));
t99 = t93 * t95;
t97 = cos(qJ(2));
t98 = t93 * t97;
t96 = cos(qJ(3));
t94 = sin(qJ(3));
t1 = [0, t101, t90 * t98 + t92 * t95, 0 (-t90 * t99 + t92 * t97) * t94 - t96 * t101, 0; 0, -t100, t90 * t95 - t92 * t98, 0 (t90 * t97 + t92 * t99) * t94 + t96 * t100, 0; 0, t93, -t91 * t97, 0, t91 * t95 * t94 - t93 * t96, 0;];
Jg_rot  = t1;
