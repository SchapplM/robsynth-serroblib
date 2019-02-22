% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:26
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRP1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:26:03
% EndTime: 2019-02-22 10:26:03
% DurationCPUTime: 0.04s
% Computational Cost: add. (59->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
t79 = qJ(1) + pkin(9);
t75 = sin(t79);
t80 = sin(qJ(5));
t87 = t75 * t80;
t81 = cos(qJ(5));
t86 = t75 * t81;
t78 = qJ(3) + pkin(10);
t76 = cos(t78);
t85 = t76 * t80;
t84 = t76 * t81;
t77 = cos(t79);
t83 = t77 * t80;
t82 = t77 * t81;
t74 = sin(t78);
t73 = t76 * t82 + t87;
t72 = -t76 * t83 + t86;
t71 = -t75 * t84 + t83;
t70 = t75 * t85 + t82;
t1 = [t71, 0, -t74 * t82, 0, t72, 0; t73, 0, -t74 * t86, 0, -t70, 0; 0, 0, t84, 0, -t74 * t80, 0; t70, 0, t74 * t83, 0, -t73, 0; t72, 0, t74 * t87, 0, t71, 0; 0, 0, -t85, 0, -t74 * t81, 0; -t75 * t74, 0, t77 * t76, 0, 0, 0; t77 * t74, 0, t75 * t76, 0, 0, 0; 0, 0, t74, 0, 0, 0;];
JR_rot  = t1;
