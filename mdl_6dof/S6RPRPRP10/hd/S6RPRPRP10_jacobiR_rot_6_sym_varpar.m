% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRP10_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_jacobiR_rot_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:48:37
% EndTime: 2019-02-26 20:48:37
% DurationCPUTime: 0.03s
% Computational Cost: add. (13->11), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
t75 = sin(qJ(3));
t76 = sin(qJ(1));
t86 = t76 * t75;
t77 = cos(qJ(5));
t85 = t76 * t77;
t74 = sin(qJ(5));
t78 = cos(qJ(3));
t84 = t78 * t74;
t83 = t78 * t77;
t79 = cos(qJ(1));
t82 = t79 * t75;
t81 = t79 * t77;
t80 = t79 * t78;
t73 = -t76 * t84 + t81;
t72 = t79 * t74 + t76 * t83;
t71 = t74 * t80 + t85;
t70 = -t76 * t74 + t77 * t80;
t1 = [-t71, 0, t74 * t86, 0, -t72, 0; t73, 0, -t74 * t82, 0, t70, 0; 0, 0, t84, 0, t75 * t77, 0; t82, 0, t76 * t78, 0, 0, 0; t86, 0, -t80, 0, 0, 0; 0, 0, -t75, 0, 0, 0; t70, 0, -t75 * t85, 0, t73, 0; t72, 0, t75 * t81, 0, t71, 0; 0, 0, -t83, 0, t75 * t74, 0;];
JR_rot  = t1;
